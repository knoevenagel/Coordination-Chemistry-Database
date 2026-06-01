"""Task registry manager for model-after task uploads (Phase 3B MVP)."""

from __future__ import annotations

import csv
import json
import re
import shutil
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Optional

TASK_ID_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{2,63}$")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _task_root(workspace_root: str | Path, task_id: str) -> Path:
    if not TASK_ID_RE.match(task_id):
        raise ValueError(f"invalid task_id: {task_id!r}")
    return Path(workspace_root).resolve() / "model_after_tasks" / task_id


def _count_rows_strict_csv(content: bytes | str, *, label: str) -> int:
    if isinstance(content, bytes):
        text = content.decode("utf-8")
    else:
        text = content
    reader = csv.DictReader(StringIO(text))
    if not (reader.fieldnames or []):
        raise ValueError(f"{label}: CSV header is required")
    rows = 0
    for row in reader:
        if any((v or "").strip() for v in row.values()):
            rows += 1
    if rows <= 0:
        raise ValueError(f"{label}: CSV must contain at least one data row")
    return rows


def create_task_from_uploads(
    *,
    task_id: str,
    workspace_root: str | Path,
    candidates_csv: bytes | str,
    positives_csv: bytes | str,
    negatives_csv: Optional[bytes | str] = None,
    metal: Optional[str] = None,
    embedding: Optional[str] = None,
    notes: str = "",
) -> Dict[str, Any]:
    task_dir = _task_root(workspace_root, task_id)
    if task_dir.exists():
        raise FileExistsError(f"task already exists: {task_dir}")

    candidate_count = _count_rows_strict_csv(candidates_csv, label="candidates")
    positive_count = _count_rows_strict_csv(positives_csv, label="positives")
    negative_count = 0
    if negatives_csv is not None and str(negatives_csv).strip() != "":
        negative_count = _count_rows_strict_csv(negatives_csv, label="negatives")

    task_dir.mkdir(parents=True, exist_ok=False)
    try:
        candidates_name = "candidates.csv"
        positives_name = "positives.csv"
        negatives_name = "negatives.csv" if negatives_csv is not None and str(negatives_csv).strip() != "" else ""

        def _write(name: str, content: bytes | str) -> None:
            p = task_dir / name
            if isinstance(content, bytes):
                p.write_bytes(content)
            else:
                p.write_text(content, encoding="utf-8")

        _write(candidates_name, candidates_csv)
        _write(positives_name, positives_csv)
        if negatives_name:
            _write(negatives_name, negatives_csv or b"")

        payload = {
            "task_id": task_id,
            "created_at": _utc_now(),
            "metal": metal or "",
            "embedding": embedding or "",
            "candidate_file": candidates_name,
            "positive_file": positives_name,
            "negative_file": negatives_name,
            "notes": notes or "",
            "candidate_count": candidate_count,
            "positive_count": positive_count,
            "negative_count": negative_count,
        }
        (task_dir / "task.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception:
        shutil.rmtree(task_dir, ignore_errors=True)
        raise

    return {
        "task_id": task_id,
        "task_dir": str(task_dir.resolve()),
        "candidate_count": candidate_count,
        "positive_count": positive_count,
        "negative_count": negative_count,
    }
