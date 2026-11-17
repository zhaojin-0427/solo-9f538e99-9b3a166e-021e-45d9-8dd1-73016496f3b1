"""SQLite 存储层测试：修订号、事务回滚（失败不改变已有判断）、重启持久化。"""

import os
import sqlite3
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.storage import Store

LOCI = [f"SNP{i}" for i in range(1, 7)]
LEFT = [("L1", [0, 0, 0, 0, 0, 0]), ("L2", [1, 1, 1, 1, 1, 1])]
RIGHT = [("R1", [0, 0, 0, 0, 0, 0]), ("R2", [1, 1, 1, 1, 1, 1])]


@pytest.fixture()
def store(tmp_path):
    return Store(str(tmp_path / "s.db"))


def test_create_and_read(store):
    ds = store.create_dataset("左", "右", LOCI, LEFT, RIGHT, optimum=0,
                              feasible=True, reason=None)
    row = store.get_dataset(ds)
    assert row["revision"] == 1 and row["feasible"] == 1 and row["optimum"] == 0
    assert store.get_samples(ds, "left") == LEFT
    assert store.get_samples(ds, "right") == RIGHT
    assert store.get_locks(ds) == {}


def test_lock_bumps_revision(store):
    ds = store.create_dataset("左", "右", LOCI, LEFT, RIGHT, 0, True, None)
    store.add_lock(ds, "L1", "R1")
    assert store.get_locks(ds) == {"L1": "R1"}
    assert store.get_dataset(ds)["revision"] == 2
    store.add_lock(ds, "L2", "R2")
    assert store.get_dataset(ds)["revision"] == 3
    store.remove_lock(ds, "L1")
    assert store.get_locks(ds) == {"L2": "R2"}
    assert store.get_dataset(ds)["revision"] == 4


def test_transaction_rollback_keeps_state(store):
    ds = store.create_dataset("左", "右", LOCI, LEFT, RIGHT, 0, True, None)
    store.add_lock(ds, "L1", "R1")  # revision -> 2
    before_rev = store.get_dataset(ds)["revision"]
    before_locks = store.get_locks(ds)

    # 模拟“校验通过后写入阶段出错”：事务中先加锁再触发约束异常
    with pytest.raises(sqlite3.Error):
        with store.transaction() as conn:
            conn.execute(
                "INSERT INTO locks (dataset_id, left_id, right_id, locked_at) "
                "VALUES (?, 'L2', 'R2', 'x')", (ds,))
            # 重复主键触发 IntegrityError
            conn.execute(
                "INSERT INTO locks (dataset_id, left_id, right_id, locked_at) "
                "VALUES (?, 'L2', 'R2', 'x')", (ds,))

    # 回滚后锁表与修订号均保持失败前状态
    assert store.get_locks(ds) == before_locks
    assert store.get_dataset(ds)["revision"] == before_rev


def test_persistence_across_store_recreate(tmp_path):
    db = str(tmp_path / "persist.db")
    s1 = Store(db)
    ds = s1.create_dataset("左", "右", LOCI, LEFT, RIGHT, 0, True, None)
    s1.add_lock(ds, "L1", "R1")
    # 重新打开同一文件（等价于进程重启）
    s2 = Store(db)
    assert s2.get_dataset(ds)["revision"] == 2
    assert s2.get_locks(ds) == {"L1": "R1"}
    assert s2.get_samples(ds, "left") == LEFT


def test_wal_mode_enabled(store):
    with store._connect() as conn:
        mode = conn.execute("PRAGMA journal_mode").fetchone()[0]
    assert mode.lower() == "wal"


def test_snapshot_reads_single_committed_state(store):
    """读快照事务期间其他连接提交新锁定，快照内读取仍看到旧的已提交状态。"""
    ds = store.create_dataset("左", "右", LOCI, LEFT, RIGHT, 0, True, None)
    with store.snapshot() as conn:
        q_locks = "SELECT COUNT(*) FROM locks WHERE dataset_id = ?"
        q_rev = "SELECT revision FROM datasets WHERE id = ?"
        before_locks = conn.execute(q_locks, (ds,)).fetchone()[0]
        before_rev = conn.execute(q_rev, (ds,)).fetchone()[0]
        # 另一连接在快照读取期间提交新锁定（WAL 下写不阻塞读）
        store.add_lock(ds, "L1", "R1")
        after_locks = conn.execute(q_locks, (ds,)).fetchone()[0]
        after_rev = conn.execute(q_rev, (ds,)).fetchone()[0]
    # 快照内两次读取一致：同一个已提交状态，看不到新提交
    assert before_locks == after_locks == 0
    assert before_rev == after_rev == 1
    # 快照结束后，新读取看到新状态
    assert store.get_locks(ds) == {"L1": "R1"}
    assert store.get_dataset(ds)["revision"] == 2


def test_load_all_bundles_returns_everything(store):
    ds1 = store.create_dataset("左", "右", LOCI, LEFT, RIGHT, 0, True, None)
    ds2 = store.create_dataset("左2", "右2", LOCI, LEFT, RIGHT, 0, True, None)
    store.add_lock(ds2, "L1", "R1")
    bundles = store.load_all_bundles()
    assert [b[0]["id"] for b in bundles] == [ds1, ds2]
    _, samples2, locks2 = bundles[1]
    assert len(samples2) == 4  # 两侧各 2 个样本
    assert {r["left_id"]: r["right_id"] for r in locks2} == {"L1": "R1"}
