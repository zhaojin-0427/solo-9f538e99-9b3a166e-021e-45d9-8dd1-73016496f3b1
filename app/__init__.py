"""Flask 应用：SNP 指纹匿名样本核对网页。

路由
----
GET  /                     任务列表
GET  /upload               上传表单
POST /upload               上传两份 CSV，校验后建任务
GET  /datasets/<id>        核对详情（候选、锁定、修订号、跨任务核查）
POST /datasets/<id>/lock   按外部证据锁定候选（带修订号乐观并发控制）
POST /datasets/<id>/unlock 解除锁定
GET  /datasets/<id>/download 全部唯一且无跨任务冲突时下载编号对照 CSV

跨任务核查
----------
详情页与下载共用一个读快照：同一院内编号在多个“已全部唯一”任务中对应的
匿名样本，原六个位点双方非空且数值不同即为矛盾；当前任务有矛盾则禁止下载。
读取失败一律不按“无冲突”处理（页面报错、下载取消）。
"""

from __future__ import annotations

import io
import json
import os
import sqlite3
from typing import Any, Dict, List, Optional

from flask import (Flask, abort, flash, redirect, render_template, request,
                   send_file, url_for)

from . import csvio
from .crosscheck import ResolvedTask, find_conflicts
from .matching import all_resolved, analyze
from .storage import Store


def create_app(db_path: Optional[str] = None) -> Flask:
    app = Flask(__name__)
    app.secret_key = os.environ.get("APP_SECRET_KEY", "dev-only-change-me")
    # 每份 CSV 最多 12 行，2MB 足够且可防止超大上传耗尽内存
    app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024
    db_path = db_path or os.environ.get(
        "APP_DB_PATH", os.path.join(os.getcwd(), "data", "snp.db")
    )
    store = Store(db_path)
    app.extensions["store"] = store

    # ----------------------------------------------------------- 工具函数
    def load_all_tasks() -> Dict[int, Dict[str, Any]]:
        """在一个读快照中加载全部任务并逐个重算分析。

        快照内所有读取来自同一个已提交状态：交错的上传/锁定提交要么整体
        可见、要么整体不可见，保证每次页面渲染或下载的“对应关系 + 冲突
        判定”基于同一份状态。读取失败抛 sqlite3.Error / ValueError，
        调用方必须按“不能确认无冲突”处理。
        """
        tasks: Dict[int, Dict[str, Any]] = {}
        for row, sample_rows, lock_rows in store.load_all_bundles():
            loci = json.loads(row["loci_json"])
            left_samples: List[Any] = []
            right_samples: List[Any] = []
            for s in sample_rows:
                item = (s["sample_id"], json.loads(s["genotype_json"]))
                (left_samples if s["side"] == "left" else right_samples
                 ).append(item)
            locks = {r["left_id"]: r["right_id"] for r in lock_rows}
            result = analyze(left_samples, right_samples, loci, locks)
            tasks[row["id"]] = {
                "row": row,
                "loci": loci,
                "left_samples": left_samples,
                "right_samples": right_samples,
                "locks": locks,
                "result": result,
            }
        return tasks

    def resolved_view(task: Dict[str, Any]) -> ResolvedTask:
        """已全部唯一任务的比对视图：院内编号 -> (匿名样本号, 指纹)。"""
        result = task["result"]
        geno_by_left = dict(task["left_samples"])
        by_right: Dict[str, Any] = {}
        for lid in result.left_ids:
            rid = result.candidates[lid][0]
            by_right[rid] = (lid, tuple(geno_by_left[lid]))
        return ResolvedTask(
            dataset_id=task["row"]["id"],
            loci=tuple(task["loci"]),
            by_right=by_right,
        )

    def cross_status(dataset_id: int,
                     tasks: Dict[int, Dict[str, Any]]) -> Dict[str, Any]:
        """当前任务的跨任务核查状态。

        status:
          not_resolved 当前任务尚未全部唯一，不参与比对；
          no_peer      没有其他已全部唯一的可比任务；
          clear        有可比任务，未发现矛盾；
          conflict     发现同号院内编号指纹矛盾。
        """
        current = tasks[dataset_id]
        base = {"peers": [], "conflicts": [], "checked": 0}
        if not all_resolved(current["result"]):
            return dict(base, status="not_resolved")
        cur_view = resolved_view(current)
        peers = [
            resolved_view(t) for tid, t in sorted(tasks.items())
            if tid != dataset_id and all_resolved(t["result"])
        ]
        if not peers:
            return dict(base, status="no_peer")
        conflicts = find_conflicts(cur_view, peers)
        checked = sum(
            1 for p in peers for rid in cur_view.by_right if rid in p.by_right
        )
        return {
            "status": "conflict" if conflicts else "clear",
            "peers": peers,
            "conflicts": conflicts,
            "checked": checked,
        }

    def read_failure(message: str):
        return render_template(
            "conflict.html",
            title="数据读取失败",
            message=message + "。读取失败不会被视为无冲突：在成功重新读取前，"
                              "相关任务的核查与下载均不可用。",
            dataset_id=None,
        ), 500

    def conflict(message: str, dataset_id: int):
        return render_template(
            "conflict.html", message=message, dataset_id=dataset_id
        ), 409

    # --------------------------------------------------------------- 首页
    @app.route("/")
    def index():
        return render_template("index.html", datasets=store.list_datasets())

    # --------------------------------------------------------------- 上传
    @app.route("/upload", methods=["GET"])
    def upload_form():
        return render_template("upload.html")

    @app.route("/upload", methods=["POST"])
    def upload_submit():
        left_file = request.files.get("left_csv")
        right_file = request.files.get("right_csv")
        left_label = (request.form.get("left_label") or "匿名侧（测序结果）").strip()
        right_label = (request.form.get("right_label") or "院内侧（库存编号）").strip()

        if left_file is None or right_file is None:
            flash("请同时选择两份 CSV 文件后再上传。", "error")
            return render_template("upload.html"), 400
        left_raw = left_file.read()
        right_raw = right_file.read()

        left_parsed, left_errors = csvio.parse_snp_csv(left_raw, left_label)
        right_parsed, right_errors = csvio.parse_snp_csv(right_raw, right_label)
        errors = left_errors + right_errors
        if errors:
            for msg in errors:
                flash(msg, "error")
            return render_template("upload.html"), 400

        left_samples, right_samples, pair_errors = csvio.validate_pair(
            left_parsed, right_parsed
        )
        if pair_errors:
            for msg in pair_errors:
                flash(msg, "error")
            return render_template("upload.html"), 400

        # 计算在数据库之外完成；任何异常都不会产生半成品数据
        result = analyze(left_samples, right_samples, left_parsed.loci)
        dataset_id = store.create_dataset(
            left_label=left_label,
            right_label=right_label,
            loci=left_parsed.loci,
            left_samples=left_samples,
            right_samples=right_samples,
            optimum=result.optimum,
            feasible=result.feasible,
            reason=result.reason,
        )
        if result.feasible:
            flash(
                "上传成功：最优不一致位点总数为 %d。请逐样本核对候选并按外部证据消歧。"
                % result.optimum,
                "info",
            )
        else:
            flash("上传成功，但该批数据无解：%s" % result.reason, "error")
        return redirect(url_for("dataset_detail", dataset_id=dataset_id))

    # --------------------------------------------------------------- 详情
    @app.route("/datasets/<int:dataset_id>")
    def dataset_detail(dataset_id: int):
        try:
            tasks = load_all_tasks()
        except (sqlite3.Error, ValueError) as exc:
            return read_failure("读取任务数据失败（%s）" % exc)
        task = tasks.get(dataset_id)
        if task is None:
            abort(404)
        result = task["result"]
        ctx = dict(task)
        ctx["all_resolved"] = all_resolved(result)
        ctx["cross"] = cross_status(dataset_id, tasks)
        # 供模板渲染的候选详情：left_id -> [{right_id, shared, mismatches, locked}]
        candidate_details: Dict[str, List[Dict[str, Any]]] = {}
        for lid in result.left_ids:
            items = []
            for rid in result.candidates.get(lid, []):
                info = result.pair_infos[(lid, rid)]
                items.append({
                    "right_id": rid,
                    "shared": info.shared,
                    "mismatches": info.mismatches,
                    "locked": task["locks"].get(lid) == rid,
                })
            candidate_details[lid] = items
        ctx["candidate_details"] = candidate_details
        return render_template("dataset.html", **ctx)

    # --------------------------------------------------------------- 锁定
    @app.route("/datasets/<int:dataset_id>/lock", methods=["POST"])
    def lock_candidate(dataset_id: int):
        left_id = (request.form.get("left_id") or "").strip()
        right_id = (request.form.get("right_id") or "").strip()
        try:
            seen_revision = int(request.form.get("revision", ""))
        except (TypeError, ValueError):
            return conflict("表单缺少修订号，请刷新页面后重试。", dataset_id)

        # 整个“校验 + 写入”放在同一个 IMMEDIATE 事务里，失败即回滚
        try:
            with store.transaction() as conn:
                row = conn.execute(
                    "SELECT * FROM datasets WHERE id = ?", (dataset_id,)
                ).fetchone()
                if row is None:
                    abort(404)
                if row["revision"] != seen_revision:
                    return conflict(
                        "页面修订号已过期（您看到的是 r%d，当前为 r%d），"
                        "可能有其他人刚提交了锁定。请刷新页面查看最新判断后再操作。"
                        % (seen_revision, row["revision"]),
                        dataset_id,
                    )

                locks = {
                    r["left_id"]: r["right_id"]
                    for r in conn.execute(
                        "SELECT left_id, right_id FROM locks WHERE dataset_id = ?",
                        (dataset_id,),
                    ).fetchall()
                }
                left_samples = [
                    (r["sample_id"], json.loads(r["genotype_json"]))
                    for r in conn.execute(
                        "SELECT sample_id, genotype_json FROM samples "
                        "WHERE dataset_id=? AND side='left' ORDER BY idx",
                        (dataset_id,),
                    ).fetchall()
                ]
                right_samples = [
                    (r["sample_id"], json.loads(r["genotype_json"]))
                    for r in conn.execute(
                        "SELECT sample_id, genotype_json FROM samples "
                        "WHERE dataset_id=? AND side='right' ORDER BY idx",
                        (dataset_id,),
                    ).fetchall()
                ]
                loci = json.loads(row["loci_json"])

                left_ids = [s[0] for s in left_samples]
                right_ids = [s[0] for s in right_samples]
                if left_id not in left_ids or right_id not in right_ids:
                    flash("锁定失败：样本号不存在于本任务中。", "error")
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))
                if left_id in locks:
                    flash(
                        "锁定失败：%s 已锁定到 %s，请先解除原锁定再改选。"
                        % (left_id, locks[left_id]),
                        "error",
                    )
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))
                used_rights = set(locks.values())
                if right_id in used_rights:
                    flash(
                        "锁定失败：%s 已被其他样本锁定，一对一对应不允许共用。"
                        % right_id,
                        "error",
                    )
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))

                trial_locks = dict(locks)
                trial_locks[left_id] = right_id
                result = analyze(left_samples, right_samples, loci, trial_locks)
                base_optimum = row["optimum"]

                if not result.feasible or result.forced_total is None:
                    flash("锁定被拒绝：%s。已有锁定保持不变。"
                          % (result.reason or "无法构成合法匹配"), "error")
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))
                if result.forced_total != base_optimum:
                    flash(
                        "锁定被拒绝：该选择会使不一致位点总数偏离原始最优值"
                        "（最优 %d，当前选择将达 %d）。已有锁定保持不变。"
                        % (base_optimum, result.forced_total),
                        "error",
                    )
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))

                conn.execute(
                    "INSERT INTO locks (dataset_id, left_id, right_id, locked_at) "
                    "VALUES (?, ?, ?, strftime('%Y-%m-%dT%H:%M:%S','now'))",
                    (dataset_id, left_id, right_id),
                )
                conn.execute(
                    "UPDATE datasets SET revision = revision + 1 WHERE id = ?",
                    (dataset_id,),
                )
        except sqlite3.Error as exc:
            flash("数据库保存失败，已有判断未被改变：%s" % exc, "error")
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))

        flash("已锁定 %s ↔ %s，其余样本候选已按最优约束重新计算。"
              % (left_id, right_id), "info")
        return redirect(url_for("dataset_detail", dataset_id=dataset_id))

    # --------------------------------------------------------------- 解锁
    @app.route("/datasets/<int:dataset_id>/unlock", methods=["POST"])
    def unlock_candidate(dataset_id: int):
        left_id = (request.form.get("left_id") or "").strip()
        try:
            seen_revision = int(request.form.get("revision", ""))
        except (TypeError, ValueError):
            return conflict("表单缺少修订号，请刷新页面后重试。", dataset_id)

        try:
            with store.transaction() as conn:
                row = conn.execute(
                    "SELECT * FROM datasets WHERE id = ?", (dataset_id,)
                ).fetchone()
                if row is None:
                    abort(404)
                if row["revision"] != seen_revision:
                    return conflict(
                        "页面修订号已过期（您看到的是 r%d，当前为 r%d），"
                        "请刷新页面后再解除锁定。"
                        % (seen_revision, row["revision"]),
                        dataset_id,
                    )
                cur = conn.execute(
                    "DELETE FROM locks WHERE dataset_id = ? AND left_id = ?",
                    (dataset_id, left_id),
                )
                if cur.rowcount == 0:
                    flash("该样本本就没有锁定，无需解除。", "error")
                    return redirect(url_for("dataset_detail",
                                            dataset_id=dataset_id))
                conn.execute(
                    "UPDATE datasets SET revision = revision + 1 WHERE id = ?",
                    (dataset_id,),
                )
        except sqlite3.Error as exc:
            flash("数据库保存失败，已有判断未被改变：%s" % exc, "error")
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))

        flash("已解除 %s 的锁定，候选已重新计算。" % left_id, "info")
        return redirect(url_for("dataset_detail", dataset_id=dataset_id))

    # --------------------------------------------------------------- 下载
    @app.route("/datasets/<int:dataset_id>/download")
    def download_mapping(dataset_id: int):
        # 对应关系与跨任务冲突判定必须来自同一个已提交状态：
        # 一次快照读完，再据此决定是否供数。
        try:
            tasks = load_all_tasks()
        except (sqlite3.Error, ValueError) as exc:
            flash("读取任务数据失败（%s），无法确认跨任务核查结果；"
                  "读取失败不按无冲突处理，本次下载已取消。" % exc, "error")
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))
        task = tasks.get(dataset_id)
        if task is None:
            abort(404)
        result = task["result"]
        if not result.feasible:
            flash("无解任务不能下载编号对照：%s" % result.reason, "error")
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))
        if not all_resolved(result):
            ambiguous = [lid for lid, c in result.candidates.items()
                         if len(c) != 1]
            flash(
                "仍有 %d 个样本存在歧义（%s），按规定不得供数；"
                "请继续按外部证据锁定消歧。"
                % (len(ambiguous), "、".join(ambiguous)),
                "error",
            )
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))

        cross = cross_status(dataset_id, tasks)
        if cross["status"] == "conflict":
            others = sorted({c.other_dataset_id
                             for c in cross["conflicts"]})
            flash(
                "跨任务核查发现 %d 处同号院内编号指纹矛盾（涉及任务 %s），"
                "按规定禁止下载编号对照；请先核实矛盾位点。"
                % (len(cross["conflicts"]),
                   "、".join("#%d" % i for i in others)),
                "error",
            )
            return redirect(url_for("dataset_detail", dataset_id=dataset_id))

        buf = io.StringIO()
        buf.write("left_sample_id,right_sample_id,mismatched_loci,shared_loci\n")
        for lid in result.left_ids:
            rid = result.candidates[lid][0]
            info = result.pair_infos[(lid, rid)]
            buf.write("%s,%s,%d,%d\n"
                      % (lid, rid, info.mismatches, info.shared))
        data = buf.getvalue().encode("utf-8-sig")
        return send_file(
            io.BytesIO(data),
            mimetype="text/csv; charset=utf-8",
            as_attachment=True,
            download_name="mapping_dataset_%d.csv" % dataset_id,
        )

    @app.after_request
    def no_store(resp):
        # 页面与下载每次都由服务器重新判定，禁止浏览器/中间缓存复用旧判断
        resp.headers["Cache-Control"] = "no-store"
        return resp

    @app.errorhandler(404)
    def not_found(_err):
        return render_template("conflict.html",
                               title="任务不存在",
                               message="任务不存在或已被删除。",
                               dataset_id=None), 404

    @app.errorhandler(413)
    def too_large(_err):
        flash("上传文件过大（单请求上限 2MB），请检查 CSV 内容。", "error")
        return render_template("upload.html"), 413

    return app
