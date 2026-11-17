"""跨任务核查：同一院内编号在多个已全部唯一任务中对应的匿名样本指纹是否矛盾。

判定规则
--------
- 只有“已全部唯一”的任务才参与比对（其匿名样本 -> 院内编号对应已确定）；
- 对同一院内编号，比较两个任务各自对应匿名样本的原六个位点指纹；
- 某位点双方均非空且数值不同 -> 记为一处矛盾；
- 任一侧缺失（空）的位点不参与判定；
- 位点按名称对齐（不同任务的位点命名/顺序可能不同，只比较同名位点）。

本模块为纯函数，不依赖 Web/数据库，便于单元测试。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class ResolvedTask:
    """一个已全部唯一任务的比对视图。"""

    dataset_id: int
    loci: Tuple[str, ...]
    # 院内编号 -> (匿名样本号, 六个位点值，与 loci 对齐)
    by_right: Dict[str, Tuple[str, Tuple[Optional[int], ...]]]


@dataclass(frozen=True)
class CrossConflict:
    """一处跨任务指纹矛盾。"""

    right_id: str          # 院内编号
    current_dataset_id: int
    other_dataset_id: int
    current_left_id: str   # 当前任务中对应的匿名样本
    other_left_id: str     # 对方任务中对应的匿名样本
    locus: str             # 位点名
    current_value: int     # 当前任务匿名样本在该位点的值
    other_value: int       # 对方任务匿名样本在该位点的值


def find_conflicts(current: ResolvedTask,
                   peers: Sequence[ResolvedTask]) -> List[CrossConflict]:
    """找出 current 与每个 peer 之间同号院内编号的指纹矛盾。

    只比较双方都已全部唯一任务中同号院内编号对应的匿名样本；
    返回按（当前任务院内编号出现顺序、对方任务号、位点顺序）排列的矛盾列表。
    """
    conflicts: List[CrossConflict] = []
    for rid, (cur_left, cur_geno) in current.by_right.items():
        for peer in sorted(peers, key=lambda p: p.dataset_id):
            hit = peer.by_right.get(rid)
            if hit is None:
                continue
            peer_left, peer_geno = hit
            peer_pos = {name: i for i, name in enumerate(peer.loci)}
            for idx, locus in enumerate(current.loci):
                j = peer_pos.get(locus)
                if j is None:
                    continue  # 对方任务没有同名位点，无法比较
                cv = cur_geno[idx] if idx < len(cur_geno) else None
                pv = peer_geno[j] if j < len(peer_geno) else None
                if cv is None or pv is None:
                    continue  # 任一侧缺失不参与判定
                if cv != pv:
                    conflicts.append(CrossConflict(
                        right_id=rid,
                        current_dataset_id=current.dataset_id,
                        other_dataset_id=peer.dataset_id,
                        current_left_id=cur_left,
                        other_left_id=peer_left,
                        locus=locus,
                        current_value=cv,
                        other_value=pv,
                    ))
    return conflicts
