"""跨任务核查测试：纯函数判定、详情页展示、下载门禁、交错提交一致性。"""

import io
import os
import sqlite3
import sys
import threading

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app import create_app
from app.crosscheck import ResolvedTask, find_conflicts

HEADER = "sample_id,SNP1,SNP2,SNP3,SNP4,SNP5,SNP6\n"
LOCI = [f"SNP{i}" for i in range(1, 7)]


def csv_bytes(lines):
    return (HEADER + "".join(lines)).encode("utf-8")


def upload(client, left_bytes, right_bytes):
    return client.post(
        "/upload",
        data={
            "left_csv": (io.BytesIO(left_bytes), "left.csv"),
            "right_csv": (io.BytesIO(right_bytes), "right.csv"),
        },
        content_type="multipart/form-data",
        follow_redirects=True,
    )


@pytest.fixture()
def app(tmp_path):
    app = create_app(db_path=str(tmp_path / "test.db"))
    app.config["TESTING"] = True
    return app


@pytest.fixture()
def client(app):
    return app.test_client()


# 典型数据 -------------------------------------------------------------------
# 任务 1：唯一，H1 ↔ L1（全 0）
T1_LEFT = csv_bytes(["L1,0,0,0,0,0,0\n"])
T1_RIGHT = csv_bytes(["H1,0,0,0,0,0,0\n"])
# 任务 2：唯一，H1 ↔ M1（全 1）—— 与任务 1 同号院内编号但指纹矛盾
T2_LEFT = csv_bytes(["M1,1,1,1,1,1,1\n"])
T2_RIGHT = csv_bytes(["H1,1,1,1,1,1,1\n"])
# 任务 3：唯一，H1 ↔ N1（全 0）—— 与任务 1 指纹一致，无矛盾
T3_LEFT = csv_bytes(["N1,0,0,0,0,0,0\n"])
T3_RIGHT = csv_bytes(["H1,0,0,0,0,0,0\n"])
# 任务 4：歧义（M1/M2 指纹相同），锁定 M1↔H1 后才全部唯一
T4_LEFT = csv_bytes(["M1,1,1,1,1,1,1\n", "M2,1,1,1,1,1,1\n"])
T4_RIGHT = csv_bytes(["H1,1,1,1,1,1,1\n", "H2,1,1,1,1,1,1\n"])


# ---------------------------------------------------------------------------
# 纯函数 find_conflicts
# ---------------------------------------------------------------------------
def rt(dataset_id, by_right, loci=None):
    return ResolvedTask(dataset_id=dataset_id,
                        loci=tuple(loci or LOCI), by_right=by_right)


def test_find_conflicts_basic_fields():
    cur = rt(1, {"H1": ("L1", (0, 0, 0, 0, 0, 0))})
    peer = rt(2, {"H1": ("M1", (1, 0, 0, 0, 0, 2))})
    conflicts = find_conflicts(cur, [peer])
    assert len(conflicts) == 2
    c = conflicts[0]
    assert (c.right_id, c.current_dataset_id, c.other_dataset_id) == ("H1", 1, 2)
    assert (c.current_left_id, c.other_left_id) == ("L1", "M1")
    assert c.locus == "SNP1" and c.current_value == 0 and c.other_value == 1
    assert conflicts[1].locus == "SNP6" and conflicts[1].other_value == 2


def test_find_conflicts_skips_missing_and_equal():
    cur = rt(1, {"H1": ("L1", (0, 0, 0, 0, 0, 0))})
    # 对方在 SNP1/SNP3/SNP5 缺失，其余相等 -> 无矛盾
    peer = rt(2, {"H1": ("M1", (None, 0, None, 0, None, 0))})
    assert find_conflicts(cur, [peer]) == []
    # 当前侧缺失同样不参与
    cur2 = rt(1, {"H1": ("L1", (None, None, None, None, None, None))})
    peer2 = rt(2, {"H1": ("M1", (1, 1, 1, 1, 1, 1))})
    assert find_conflicts(cur2, [peer2]) == []


def test_find_conflicts_aligns_loci_by_name():
    # 对方位点顺序不同且少一个同名位点：只按同名位点比较
    cur = rt(1, {"H1": ("L1", (0, 1, 0, 0, 0, 0))})
    peer_loci = ["SNP2", "SNP1", "SNP3", "SNP4", "SNP5", "OTHER"]
    peer = rt(2, {"H1": ("M1", (0, 1, 0, 0, 0, 9))}, loci=peer_loci)
    conflicts = find_conflicts(cur, [peer])
    # SNP1: 0 vs 1 矛盾；SNP2: 1 vs 0 矛盾；SNP6 对方无同名位点不参与
    assert [(c.locus, c.current_value, c.other_value) for c in conflicts] == [
        ("SNP1", 0, 1), ("SNP2", 1, 0),
    ]


def test_find_conflicts_multiple_peers_and_unshared_right():
    cur = rt(1, {"H1": ("L1", (0, 0, 0, 0, 0, 0)),
                 "H9": ("L9", (0, 0, 0, 0, 0, 0))})
    peer_a = rt(2, {"H1": ("M1", (1, 1, 1, 1, 1, 1))})
    peer_b = rt(3, {"H1": ("N1", (0, 0, 0, 0, 0, 0)),
                    "H9": ("N9", (2, 2, 2, 2, 2, 2))})
    peer_c = rt(4, {"H7": ("X1", (1, 1, 1, 1, 1, 1))})  # 无同号院内编号
    conflicts = find_conflicts(cur, [peer_a, peer_b, peer_c])
    assert len(conflicts) == 12  # H1 vs 任务2 六处 + H9 vs 任务3 六处
    assert {c.other_dataset_id for c in conflicts} == {2, 3}
    # 任务 2 的矛盾排在任务 3 之前（按任务号排序）
    assert [c.other_dataset_id for c in conflicts] == [2] * 6 + [3] * 6


# ---------------------------------------------------------------------------
# 详情页展示与下载门禁
# ---------------------------------------------------------------------------
def test_conflict_listed_and_download_blocked(client):
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, T2_LEFT, T2_RIGHT)
    for did in (1, 2):
        body = client.get("/datasets/%d" % did).get_data(as_text=True)
        assert "跨任务核查" in body
        assert "发现 6 处矛盾" in body
        assert "L1" in body and "M1" in body and "SNP1" in body
        assert "#1" in body and "#2" in body  # 列出两任务
        assert "下载（存在跨任务矛盾）" in body
        resp = client.get("/datasets/%d/download" % did,
                          follow_redirects=True)
        assert "禁止下载编号对照" in resp.get_data(as_text=True)
        assert resp.status_code == 200  # 重定向回详情页，不是 CSV


def test_same_fingerprint_no_conflict_allows_download(client):
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, T3_LEFT, T3_RIGHT)
    body = client.get("/datasets/1").get_data(as_text=True)
    assert "未发现矛盾" in body
    assert "#2" in body  # 提示已与任务 2 比对
    resp = client.get("/datasets/1/download")
    assert resp.status_code == 200
    assert "L1,H1,0,6" in resp.data.decode("utf-8-sig")


def test_no_peer_task_prompt(client):
    upload(client, T1_LEFT, T1_RIGHT)
    body = client.get("/datasets/1").get_data(as_text=True)
    assert "没有可比任务" in body
    # 没有可比任务 -> 沿用原规则，可下载
    assert client.get("/datasets/1/download").status_code == 200


def test_not_resolved_task_not_comparable(client):
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, T4_LEFT, T4_RIGHT)  # 歧义任务
    body = client.get("/datasets/2").get_data(as_text=True)
    assert "暂不参与" in body
    # 任务 1 视角：没有已全部唯一的可比任务
    body1 = client.get("/datasets/1").get_data(as_text=True)
    assert "没有可比任务" in body1
    assert client.get("/datasets/1/download").status_code == 200


def test_missing_values_not_counted_as_conflict(client):
    # 任务 5：H1 ↔ M1，M1 只在 SNP1/SNP3/SNP5 非空且与 L1 相同 -> 无矛盾
    left = csv_bytes(["M1,0,,0,,0,\n"])
    right = csv_bytes(["H1,0,0,0,0,0,0\n"])
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, left, right)
    body = client.get("/datasets/1").get_data(as_text=True)
    assert "未发现矛盾" in body
    assert client.get("/datasets/1/download").status_code == 200


def test_missing_counterpart_allows_partial_conflict(client):
    # M1 在 SNP1 与 L1 不同、SNP2/4/6 缺失 -> 仅 SNP1 一处矛盾
    left = csv_bytes(["M1,1,,0,,0,\n"])
    right = csv_bytes(["H1,0,0,0,0,0,0\n"])  # M1 与 H1：共 3 位点、异 1，仍可配
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, left, right)
    body = client.get("/datasets/1").get_data(as_text=True)
    assert "发现 1 处矛盾" in body
    assert "SNP1" in body


# ---------------------------------------------------------------------------
# 交错提交：上传/锁定后才开始的下载必须反映新冲突
# ---------------------------------------------------------------------------
def test_download_after_upload_reflects_new_conflict(client):
    upload(client, T1_LEFT, T1_RIGHT)
    # 上传任务 2 之前：任务 1 可下载
    assert client.get("/datasets/1/download").status_code == 200
    # 上传提交之后才开始的下载必须反映新冲突
    upload(client, T2_LEFT, T2_RIGHT)
    resp = client.get("/datasets/1/download")
    assert resp.status_code == 302  # 被拦截重定向，不再是 CSV
    resp2 = client.get("/datasets/2/download")
    assert resp2.status_code == 302


def test_download_after_lock_and_unlock_reflects_conflict(client):
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, T4_LEFT, T4_RIGHT)  # 歧义：任务 1 仍可下载
    assert client.get("/datasets/1/download").status_code == 200
    # 锁定使任务 2 全部唯一（M1↔H1，M2↔H2），与任务 1 冲突
    client.post("/datasets/2/lock",
                data={"revision": 1, "left_id": "M1", "right_id": "H1"})
    assert client.get("/datasets/1/download").status_code == 302
    assert client.get("/datasets/2/download").status_code == 302
    # 解锁后任务 2 不再是可比任务，任务 1 恢复可下载
    client.post("/datasets/2/unlock",
                data={"revision": 2, "left_id": "M1"})
    assert client.get("/datasets/1/download").status_code == 200


def test_locks_and_fingerprints_unchanged_by_crosscheck(client, app):
    """跨任务核查只读：原始指纹与锁定不被改动。"""
    upload(client, T1_LEFT, T1_RIGHT)
    upload(client, T4_LEFT, T4_RIGHT)
    client.post("/datasets/2/lock",
                data={"revision": 1, "left_id": "M1", "right_id": "H1"})
    client.get("/datasets/1")            # 触发核查
    client.get("/datasets/1/download")   # 触发下载门禁
    store = app.extensions["store"]
    assert store.get_locks(2) == {"M1": "H1"}
    assert store.get_samples(1, "left") == [("L1", [0, 0, 0, 0, 0, 0])]
    assert store.get_samples(2, "left")[0] == ("M1", [1, 1, 1, 1, 1, 1])


# ---------------------------------------------------------------------------
# 读取失败不能当作无冲突
# ---------------------------------------------------------------------------
def test_read_failure_fails_closed(client, app, monkeypatch):
    upload(client, T1_LEFT, T1_RIGHT)
    store = app.extensions["store"]

    def broken_snapshot():
        raise sqlite3.OperationalError("disk I/O error")

    monkeypatch.setattr(store, "snapshot", broken_snapshot)
    # 页面：报读取失败，绝不显示“未发现矛盾”
    resp = client.get("/datasets/1")
    assert resp.status_code == 500
    body = resp.get_data(as_text=True)
    assert "读取失败" in body and "未发现矛盾" not in body
    # 下载：取消而不是放出 CSV
    resp = client.get("/datasets/1/download")
    assert resp.status_code == 302
    assert "attachment" not in resp.headers.get("Content-Disposition", "")


def test_corrupt_fingerprint_fails_closed(client, app):
    upload(client, T1_LEFT, T1_RIGHT)
    store = app.extensions["store"]
    # 直接破坏库中指纹数据，模拟读取到的数据无法解析
    with store.transaction() as conn:
        conn.execute("UPDATE samples SET genotype_json = 'not-json'")
    resp = client.get("/datasets/1")
    assert resp.status_code == 500
    assert "未发现矛盾" not in resp.get_data(as_text=True)
    resp = client.get("/datasets/1/download")
    assert resp.status_code == 302


# ---------------------------------------------------------------------------
# 交错下载：每次下载的对应与冲突判定来自同一已提交状态
# ---------------------------------------------------------------------------
def test_downloads_consistent_under_concurrent_locks(app):
    client = app.test_client()
    ambig_left = csv_bytes([
        "L1,0,0,0,0,0,0\n", "L2,0,0,0,0,0,0\n", "L3,2,2,2,2,2,2\n"])
    ambig_right = csv_bytes([
        "R1,0,0,0,0,0,0\n", "R2,0,0,0,0,0,0\n", "R3,2,2,2,2,2,2\n"])
    upload(client, ambig_left, ambig_right)
    store = app.extensions["store"]
    stop = threading.Event()
    failures = []

    def writer():
        # 另一“提交者”不断锁定/解锁（直接走存储层提交事务）
        while not stop.is_set():
            try:
                store.add_lock(1, "L1", "R2")
            except sqlite3.Error:
                pass
            try:
                store.remove_lock(1, "L1")
            except sqlite3.Error:
                pass

    t = threading.Thread(target=writer)
    t.start()
    try:
        dl = app.test_client()
        for _ in range(40):
            resp = dl.get("/datasets/1/download")
            try:
                assert resp.status_code in (200, 302)
                if resp.status_code == 200:
                    text = resp.data.decode("utf-8-sig")
                    rows = [ln.split(",")
                            for ln in text.strip().splitlines()[1:]]
                    # 任何一次成功的下载都必须是同一已提交状态下的
                    # 完整一对一对应：L1 锁定到 R2 后的唯一解
                    assert [r[0] for r in rows] == ["L1", "L2", "L3"]
                    assert {(r[0], r[1]) for r in rows} == {
                        ("L1", "R2"), ("L2", "R1"), ("L3", "R3")}
            except AssertionError as exc:  # pragma: no cover
                failures.append(str(exc))
    finally:
        stop.set()
        t.join()
    assert not failures
