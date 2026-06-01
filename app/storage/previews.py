"""Read-only preview helpers for Phase 3A Streamlit pages."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import sqlite3

from . import repositories


def _empty_df() -> pd.DataFrame:
    return pd.DataFrame()


def read_csv_preview(path: str | Path, limit: int = 100) -> Dict[str, Any]:
    p = Path(path).expanduser()
    if not p.is_file():
        return {
            "ok": False,
            "path": str(p),
            "data": _empty_df(),
            "error": f"file not found: {p}",
        }
    try:
        df = pd.read_csv(p, nrows=max(1, int(limit)))
        return {"ok": True, "path": str(p), "data": df, "error": None}
    except Exception as exc:
        return {
            "ok": False,
            "path": str(p),
            "data": _empty_df(),
            "error": str(exc),
        }


def read_json_preview(path: str | Path) -> Dict[str, Any]:
    p = Path(path).expanduser()
    if not p.is_file():
        return {
            "ok": False,
            "path": str(p),
            "data": {},
            "error": f"file not found: {p}",
        }
    try:
        with open(p, "r", encoding="utf-8") as f:
            data = json.load(f)
        return {"ok": True, "path": str(p), "data": data, "error": None}
    except Exception as exc:
        return {"ok": False, "path": str(p), "data": {}, "error": str(exc)}


def resolve_task_file(task_row: Dict[str, Any], which: str) -> Optional[Path]:
    task_dir = task_row.get("task_dir")
    if not task_dir:
        return None
    base = Path(task_dir)
    mapping = {
        "candidates": task_row.get("candidate_file"),
        "positives": task_row.get("positive_file"),
        "negatives": task_row.get("negative_file"),
    }
    if which not in mapping:
        raise ValueError(f"unsupported task file type: {which}")
    rel = mapping[which]
    if not rel:
        return None
    return base / str(rel)


def preview_ga_version(
    conn: sqlite3.Connection,
    ga_set_id: str,
    ga_version_id: str,
    limit: int = 100,
) -> Dict[str, Any]:
    row = repositories.get_ga_version(conn, ga_set_id, ga_version_id)
    if not row:
        return {"ok": False, "error": f"ga version not found: {ga_set_id}/{ga_version_id}", "data": _empty_df()}
    ga_csv_path = row.get("ga_csv_path")
    if not ga_csv_path:
        return {"ok": False, "error": "ga_csv_path is empty", "data": _empty_df()}
    out = read_csv_preview(ga_csv_path, limit=limit)
    out["ga_version"] = row
    return out


def preview_task_files(task_row: Dict[str, Any], limit: int = 100) -> Dict[str, Dict[str, Any]]:
    result: Dict[str, Dict[str, Any]] = {}
    for which in ("candidates", "positives", "negatives"):
        path = resolve_task_file(task_row, which)
        if path is None:
            result[which] = {
                "ok": False,
                "path": None,
                "data": _empty_df(),
                "error": f"{which} file not configured",
            }
            continue
        result[which] = read_csv_preview(path, limit=limit)
    return result


def preview_model_result(result_row: Dict[str, Any], limit: int = 100) -> Dict[str, Dict[str, Any]]:
    return {
        "ranking": read_csv_preview(result_row.get("ranking_path") or "", limit=limit),
        "metrics": read_json_preview(result_row.get("metrics_path") or ""),
        "report": read_json_preview(result_row.get("report_path") or ""),
    }
