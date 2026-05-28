"""Per-run training index/config preparation (Phase 1C)."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List

from .run_context import RunContext

WEB_ROOT = Path(__file__).resolve().parents[2]
INDEX_TEMPLATE = WEB_ROOT / "app" / "templates" / "training_index_template.json"
CONFIG_TEMPLATE = WEB_ROOT / "app" / "templates" / "training_config_template.yaml"

def default_training_device() -> str:
    """Training always uses GPU; override with CHEMDB_TRAINING_DEVICE (e.g. cuda:0)."""
    return os.environ.get("CHEMDB_TRAINING_DEVICE", "cuda").strip() or "cuda"


PREPARE_INDEX_INPUTS = [
    "tmp/step13_kl_nl_samples.csv",
    "tmp/m_l3_pairs.csv",
    "tmp/l3_gac.json",
    "data/L3_embedding/L3_embedding_ecfp.npz",
    "data/metal_embedding/element_features_zscore.csv",
]


def _abs_under_run(ctx: RunContext, rel: str) -> Path:
    p = (ctx.run_root / rel).resolve()
    if not ctx.resolve_under_run(p):
        raise ValueError(f"Path escapes run_root: {p}")
    return p


def _check_files(ctx: RunContext, rels: List[str]) -> Dict[str, Any]:
    missing: List[str] = []
    checked: List[str] = []
    for rel in rels:
        path = _abs_under_run(ctx, rel)
        checked.append(str(path))
        if not path.is_file():
            missing.append(rel)
    return {"ok": len(missing) == 0, "missing": missing, "checked": checked}


def prepare_training_index(ctx: RunContext) -> Dict[str, Any]:
    """Write training/index.json with absolute paths under run_root."""
    chk = _check_files(ctx, PREPARE_INDEX_INPUTS)
    if not chk["ok"]:
        return {"ok": False, "error": f"missing inputs: {chk['missing']}", "input_check": chk}

    run_root = str(ctx.run_root.resolve())
    samples = str(_abs_under_run(ctx, "tmp/step13_kl_nl_samples.csv"))
    m_l3 = str(_abs_under_run(ctx, "tmp/m_l3_pairs.csv"))
    l3_gac = str(_abs_under_run(ctx, "tmp/l3_gac.json"))
    lig_path = str(_abs_under_run(ctx, "data/L3_embedding/L3_embedding_ecfp.npz"))
    lig_dir = str(_abs_under_run(ctx, "data/L3_embedding"))
    metal_path = str(_abs_under_run(ctx, "data/metal_embedding/element_features_zscore.csv"))

    template = INDEX_TEMPLATE.read_text(encoding="utf-8")
    rendered = (
        template.replace("{run_root}", run_root)
        .replace("{samples_csv}", samples)
        .replace("{m_l3_pairs}", m_l3)
        .replace("{l3_gac}", l3_gac)
        .replace("{ligand_embedding_path}", lig_path)
        .replace("{ligand_embedding_dir}", lig_dir)
        .replace("{metal_embedding_path}", metal_path)
    )
    index = json.loads(rendered)
    out_path = ctx.training_dir / "index.json"
    ctx.training_dir.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)

    return {
        "ok": True,
        "index_path": str(out_path.resolve()),
        "ligand_embedding": index.get("ligand_embedding"),
    }


def prepare_training_config(ctx: RunContext, *, integration: bool = False) -> Dict[str, Any]:
    """Write training/config.yaml with absolute paths under run_root."""
    index_path = ctx.training_dir / "index.json"
    if not index_path.is_file():
        return {"ok": False, "error": f"missing {index_path}"}

    training_dir = str(ctx.training_dir.resolve())
    index_abs = str(index_path.resolve())
    ckpt_dir = str((ctx.training_dir / "ckpts").resolve())
    device = default_training_device()
    if integration:
        # Smoke / integration: single hyperparameter set, run through epoch 3 (no sweep).
        epochs = int(os.environ.get("CHEMDB_TRAINING_INTEGRATION_EPOCHS", "3"))
        batch_size = 4
        early_stop_patience = 0
        scheduler = "null"
    else:
        epochs = 50
        batch_size = 32
        early_stop_patience = 8
        scheduler = "reduce_on_plateau"

    template = CONFIG_TEMPLATE.read_text(encoding="utf-8")
    rendered = (
        template.replace("{training_dir}", training_dir)
        .replace("{index_path}", index_abs)
        .replace("{ckpt_dir}", ckpt_dir)
        .replace("{epochs}", str(epochs))
        .replace("{batch_size}", str(batch_size))
        .replace("{device}", device)
        .replace("{early_stop_patience}", str(early_stop_patience))
        .replace("{scheduler}", scheduler)
    )
    out_path = ctx.training_dir / "config.yaml"
    out_path.write_text(rendered, encoding="utf-8")
    return {
        "ok": True,
        "config_path": str(out_path.resolve()),
        "integration": integration,
    }
