"""Load RankModel and embeddings from a Phase 1C model run."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import torch
import yaml

from training import model as model_module
from training.data import load_ligand_embedding, load_metal_embedding


@dataclass(frozen=True)
class RunModelBundle:
    model_run_root: Path
    model_id: str
    index_path: Path
    config_path: Path
    ckpt_path: Path
    index: Dict[str, Any]
    config: Dict[str, Any]
    device: torch.device
    model: torch.nn.Module
    ligand_lookup: Dict[str, Any]
    d_l: int
    metal_lookup: Dict[str, Any]
    d_m: int
    embedding_backend: str


def _resolve_under_root(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def resolve_device(device: str = "cpu") -> torch.device:
    name = (device or "cpu").strip().lower()
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if name.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(f"device={device} but torch.cuda.is_available() is False")
    return torch.device(device)


def load_run_bundle(
    model_run_root: str | Path,
    *,
    checkpoint: str | Path | None = None,
    device: str = "cpu",
    model_id: str | None = None,
) -> RunModelBundle:
    root = Path(model_run_root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"model_run_root not found: {root}")

    index_path = root / "training" / "index.json"
    config_path = root / "training" / "config.yaml"
    if not index_path.is_file():
        raise FileNotFoundError(f"Missing index.json: {index_path}")
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing config.yaml: {config_path}")

    with open(index_path, "r", encoding="utf-8") as f:
        index = json.load(f)
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    le = index.get("ligand_embedding") or {}
    backend = (le.get("backend") or config.get("data", {}).get("embedding") or "").strip().lower()
    if backend != "ecfp":
        raise ValueError(f"Phase 1D only supports ecfp, got backend={backend!r}")

    lig_path = Path(le["path"]).resolve()
    metal_path = Path((index.get("metal_embedding") or {})["path"]).resolve()
    for p, label in ((lig_path, "ligand npz"), (metal_path, "metal csv")):
        if not p.is_file():
            raise FileNotFoundError(f"Missing {label}: {p}")
        if not _resolve_under_root(p, root):
            raise ValueError(f"{label} must be under model_run_root: {p}")

    if checkpoint is None:
        ckpt_path = root / "training" / "ckpts" / "best.pt"
    else:
        ckpt_path = Path(checkpoint).resolve()
        if not ckpt_path.is_file():
            raise FileNotFoundError(f"checkpoint not found: {ckpt_path}")
        if not _resolve_under_root(ckpt_path, root):
            raise ValueError(f"checkpoint must be under model_run_root: {ckpt_path}")

    dev = resolve_device(device)
    ligand_lookup, d_l = load_ligand_embedding(lig_path)
    metal_lookup, d_m = load_metal_embedding(metal_path)

    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    m_cfg = (ckpt.get("config") or config).get("model") or {}
    state = ckpt.get("model_state_dict", ckpt)
    rank_model = model_module.RankModel(
        d_m=d_m,
        d_l=d_l,
        activation=m_cfg.get("activation", "gelu"),
        dropout=float(m_cfg.get("dropout", 0.1)),
    ).to(dev)
    rank_model.load_state_dict(state)
    rank_model.eval()

    mid = model_id or root.name
    return RunModelBundle(
        model_run_root=root,
        model_id=mid,
        index_path=index_path,
        config_path=config_path,
        ckpt_path=ckpt_path,
        index=index,
        config=config,
        device=dev,
        model=rank_model,
        ligand_lookup=ligand_lookup,
        d_l=d_l,
        metal_lookup=metal_lookup,
        d_m=d_m,
        embedding_backend=backend,
    )
