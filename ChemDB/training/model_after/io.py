"""Task input I/O: task.json + CSV schemas."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class TaskRow:
    row_id: str
    smiles: str
    did: str
    name: str
    source: str

    def has_embedding_key(self) -> bool:
        return bool(self.smiles.strip()) or bool(self.did.strip())


@dataclass
class TaskBundle:
    task_dir: Path
    task_id: str
    metal: str
    embedding: str
    candidates: List[TaskRow]
    positives: List[TaskRow]
    negatives: List[TaskRow] = field(default_factory=list)
    notes: str = ""
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def positive_ids(self) -> set[str]:
        return {r.row_id for r in self.positives}

    @property
    def negative_ids(self) -> set[str]:
        return {r.row_id for r in self.negatives}


def _read_csv_rows(path: Path, id_col: str) -> List[TaskRow]:
    rows: List[TaskRow] = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return rows
        for i, raw in enumerate(reader):
            rid = (raw.get(id_col) or raw.get(f"{id_col.split('_')[0]}_id") or "").strip()
            if not rid:
                rid = str(i)
            smiles = (raw.get("smiles") or raw.get("SMILES") or "").strip()
            did = (raw.get("did") or raw.get("DID") or "").strip()
            name = (raw.get("name") or "").strip()
            source = (raw.get("source") or "").strip()
            if not smiles and not did:
                continue
            rows.append(TaskRow(row_id=rid, smiles=smiles, did=did, name=name, source=source))
    return rows


def load_task_bundle(task_dir: str | Path) -> TaskBundle:
    root = Path(task_dir).resolve()
    task_json = root / "task.json"
    if not task_json.is_file():
        raise FileNotFoundError(f"task.json not found: {task_json}")

    with open(task_json, "r", encoding="utf-8") as f:
        spec = json.load(f)

    task_id = str(spec.get("task_id") or root.name)
    metal = str(spec.get("metal") or "").strip()
    if not metal:
        raise ValueError("task.json must specify metal")
    embedding = str(spec.get("embedding") or "ecfp").strip().lower()
    if embedding != "ecfp":
        raise ValueError(f"Only ecfp supported, got {embedding}")

    cand_file = root / str(spec.get("candidate_file") or "candidates.csv")
    pos_file = root / str(spec.get("positive_file") or "positives.csv")
    neg_file = root / str(spec.get("negative_file") or "negatives.csv")

    if not cand_file.is_file():
        raise FileNotFoundError(f"candidates file not found: {cand_file}")
    if not pos_file.is_file():
        raise FileNotFoundError(f"positives file not found: {pos_file}")

    candidates = _read_csv_rows(cand_file, "candidate_id")
    positives = _read_csv_rows(pos_file, "positive_id")
    negatives: List[TaskRow] = []
    if neg_file.is_file():
        negatives = _read_csv_rows(neg_file, "negative_id")

    if not candidates:
        raise ValueError("candidates.csv has no valid rows (need smiles or did)")
    if not positives:
        raise ValueError("positives.csv has no valid rows (need smiles or did)")

    return TaskBundle(
        task_dir=root,
        task_id=task_id,
        metal=metal,
        embedding=embedding,
        candidates=candidates,
        positives=positives,
        negatives=negatives,
        notes=str(spec.get("notes") or ""),
        raw=spec,
    )


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
