"""Workspace-level GA sets, run binding, and apply-ga tooling (Phase 1E)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .run_context import RunContext

WEB_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORKSPACE_ROOT = WEB_ROOT / "workspace"
WORKSPACE_ROOT = DEFAULT_WORKSPACE_ROOT
GA_SETS_ROOT = WORKSPACE_ROOT / "ga_sets"
CHEMDB_ROOT = WEB_ROOT / "ChemDB"
STEP4_5_SCRIPT = CHEMDB_ROOT / "src" / "step4_5.py"

GA_BINDING_MANIFEST = "manifests/ga_binding.json"
GA_STALE_MANIFEST = "manifests/ga_stale_report.json"
MATERIALIZED_GA_REL = "tmp/GA_with_id.csv"

GA_SET_ID_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{2,63}$")
VERSION_ID_RE = re.compile(r"^v\d{3}$")

STALE_DOWNSTREAM_REL: List[str] = [
    "tmp/ligand_with_gac.csv",
    "tmp/IRL_filtered.csv",
    "tmp/IRL_filtered_cleaned.csv",
    "tmp/step4_5_stats.json",
    "tmp/fragments.csv",
    "tmp/step6_7_stats.json",
    "tmp/l3_l5.json",
    "tmp/l5_l3.json",
    "tmp/l5_freq_weight.json",
    "tmp/l3_gac.json",
    "tmp/m_l3_pairs.csv",
    "tmp/metal_l3_index.csv",
    "tmp/step12_stats.json",
    "tmp/step13_kl_nl_samples.csv",
    "tmp/step13_kl_nl_samples.stats.csv",
    "tmp/step13_stats.json",
    "data/L3_embedding/L3_embedding_ecfp.npz",
    "training/index.json",
    "training/config.yaml",
    "training/split_index.json",
    "training/train_records.pkl",
    "training/val_records.pkl",
    "training/test_records.pkl",
    "training/ckpts/best.pt",
    "training/ckpts/last.pt",
    "training/ckpts/history.json",
]

STALE_NEO4J_GLOB = "tmp/neo4j_*.csv"
STALE_TRAINING_GLOB = "training/ckpts/*.pt"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def _read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def resolve_workspace_root(workspace_root: str | Path | None = None) -> Path:
    if workspace_root is not None:
        return Path(workspace_root).expanduser().resolve()
    return WORKSPACE_ROOT.resolve()


def ga_sets_root(workspace_root: str | Path | None = None) -> Path:
    return resolve_workspace_root(workspace_root) / "ga_sets"


def ga_set_dir(ga_set_id: str, *, workspace_root: str | Path | None = None) -> Path:
    if not GA_SET_ID_RE.match(ga_set_id):
        raise ValueError(f"invalid ga_set_id: {ga_set_id!r}")
    return ga_sets_root(workspace_root) / ga_set_id


def version_dir(
    ga_set_id: str,
    ga_version_id: str,
    *,
    workspace_root: str | Path | None = None,
) -> Path:
    if not VERSION_ID_RE.match(ga_version_id):
        raise ValueError(f"invalid ga_version_id: {ga_version_id!r} (expected v001, v002, ...)")
    return ga_set_dir(ga_set_id, workspace_root=workspace_root) / "versions" / ga_version_id


def _canonicalize_smiles(smiles: str) -> str:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError(f"invalid GA_SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol, canonical=True)


def _generate_ga_id(smiles: str) -> str:
    sys.path.insert(0, str(CHEMDB_ROOT / "src"))
    from step4_5 import generate_ga_id

    gid = generate_ga_id(smiles)
    if not gid:
        raise ValueError(f"could not generate GA_ID for {smiles!r}")
    return gid


def ga_csv_checksum(csv_path: Path) -> Tuple[str, int]:
    """SHA256 of canonical two-column GA table (sorted by GA_SMILES)."""
    rows: List[Dict[str, str]] = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "GA_SMILES" not in (reader.fieldnames or []) or "GA_ID" not in (reader.fieldnames or []):
            raise ValueError(f"{csv_path}: requires GA_SMILES and GA_ID columns")
        for row in reader:
            smi = (row.get("GA_SMILES") or "").strip()
            gid = (row.get("GA_ID") or "").strip()
            if smi and gid:
                rows.append({"GA_SMILES": smi, "GA_ID": gid})
    rows.sort(key=lambda r: r["GA_SMILES"])
    buf = StringIO()
    w = csv.DictWriter(buf, fieldnames=["GA_SMILES", "GA_ID"], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    digest = hashlib.sha256(buf.getvalue().encode("utf-8")).hexdigest()
    return f"sha256:{digest}", len(rows)


def materialize_ga_csv(source: Path, dest: Path) -> None:
    """Copy GA version to run snapshot (GA_SMILES, GA_ID only)."""
    rows: List[Dict[str, str]] = []
    with open(source, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smi = (row.get("GA_SMILES") or "").strip()
            gid = (row.get("GA_ID") or "").strip()
            if not smi:
                continue
            if not gid:
                gid = _generate_ga_id(smi)
            rows.append({"GA_SMILES": smi, "GA_ID": gid})
    if not rows:
        raise ValueError(f"no GA rows in {source}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["GA_SMILES", "GA_ID"])
        w.writeheader()
        w.writerows(rows)


def _next_version_id(ga_set_id: str, *, workspace_root: str | Path | None = None) -> str:
    versions_root = ga_set_dir(ga_set_id, workspace_root=workspace_root) / "versions"
    if not versions_root.is_dir():
        return "v001"
    nums = []
    for p in versions_root.iterdir():
        if p.is_dir() and VERSION_ID_RE.match(p.name):
            nums.append(int(p.name[1:]))
    n = max(nums, default=0) + 1
    return f"v{n:03d}"


def _list_version_ids(ga_set_id: str, *, workspace_root: str | Path | None = None) -> List[str]:
    root = ga_set_dir(ga_set_id, workspace_root=workspace_root) / "versions"
    if not root.is_dir():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir() and VERSION_ID_RE.match(p.name))


def _run_step4_5(
    *,
    mode: str,
    input_dir: Path,
    output_dir: Path,
    workers: Optional[int] = None,
    limit: int = -1,
) -> None:
    cmd = [
        sys.executable,
        str(STEP4_5_SCRIPT),
        "--mode",
        mode,
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
    ]
    if workers is not None:
        cmd.extend(["--workers", str(workers)])
    if limit > 0:
        cmd.extend(["--limit", str(limit)])
    env = dict(__import__("os").environ)
    parts = [str(CHEMDB_ROOT), str(CHEMDB_ROOT / "src")]
    if env.get("PYTHONPATH"):
        parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = __import__("os").pathsep.join(parts)
    proc = subprocess.run(cmd, cwd=str(CHEMDB_ROOT), env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"step4_5 --mode {mode} failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        )


def _normalize_ga_rows(rows: List[Dict[str, str]], *, generate_missing_id: bool) -> List[Dict[str, str]]:
    seen_smiles: Dict[str, str] = {}
    seen_ids: Dict[str, str] = {}
    out: List[Dict[str, str]] = []
    for row in rows:
        raw_smi = (row.get("GA_SMILES") or "").strip()
        if not raw_smi:
            continue
        smi = _canonicalize_smiles(raw_smi)
        gid = (row.get("GA_ID") or "").strip()
        if not gid:
            if not generate_missing_id:
                raise ValueError(f"missing GA_ID for GA_SMILES={smi!r}")
            gid = _generate_ga_id(smi)
        if smi in seen_smiles:
            continue
        if gid in seen_ids and seen_ids[gid] != smi:
            raise ValueError(f"duplicate GA_ID {gid!r} for different SMILES")
        seen_smiles[smi] = gid
        seen_ids[gid] = smi
        out.append({"GA_SMILES": smi, "GA_ID": gid})
    if not out:
        raise ValueError("no valid GA rows after normalization")
    return out


def _write_ga_version_files(
    ga_set_id: str,
    ga_version_id: str,
    rows: List[Dict[str, str]],
    *,
    parent_version_id: Optional[str],
    source: Dict[str, Any],
    notes: str = "",
    workspace_root: str | Path | None = None,
) -> Path:
    vdir = version_dir(ga_set_id, ga_version_id, workspace_root=workspace_root)
    if vdir.exists():
        raise FileExistsError(f"version already exists: {vdir}")
    vdir.mkdir(parents=True, exist_ok=False)
    csv_path = vdir / "GA_with_id.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["GA_SMILES", "GA_ID"])
        w.writeheader()
        w.writerows(rows)
    checksum, num_ga = ga_csv_checksum(csv_path)
    _write_json(
        vdir / "ga_version.json",
        {
            "ga_set_id": ga_set_id,
            "ga_version_id": ga_version_id,
            "parent_version_id": parent_version_id,
            "created_at": _utc_now(),
            "created_by": "ga_registry_manager",
            "source": source,
            "num_ga": num_ga,
            "checksum": checksum,
            "notes": notes,
            "columns": ["GA_SMILES", "GA_ID"],
        },
    )
    return csv_path


def _ensure_ga_set_json(
    ga_set_id: str,
    *,
    source: Dict[str, Any],
    name: str = "",
    workspace_root: str | Path | None = None,
) -> None:
    sdir = ga_set_dir(ga_set_id, workspace_root=workspace_root)
    sdir.mkdir(parents=True, exist_ok=True)
    meta_path = sdir / "ga_set.json"
    now = _utc_now()
    if meta_path.is_file():
        meta = _read_json(meta_path)
        meta["updated_at"] = now
        if source:
            meta["source"] = source
        _write_json(meta_path, meta)
    else:
        _write_json(
            meta_path,
            {
                "ga_set_id": ga_set_id,
                "name": name or ga_set_id,
                "description": "",
                "created_at": now,
                "updated_at": now,
                "tags": [],
                "source": source,
            },
        )


def _collect_stale_files(
    ctx: RunContext,
    *,
    workspace_root: str | Path | None = None,
) -> List[Dict[str, Any]]:
    stale: List[Dict[str, Any]] = []
    for rel in STALE_DOWNSTREAM_REL:
        path = ctx.run_root / rel
        if path.is_file():
            stale.append({"path": rel, "exists": True, "reason": "downstream_of_ga"})
    for path in sorted(ctx.tmp_dir.glob("neo4j_*.csv")):
        rel = str(path.relative_to(ctx.run_root))
        stale.append({"path": rel, "exists": True, "reason": "downstream_of_ga"})
    for path in sorted((ctx.training_dir / "ckpts").glob("*.pt")) if (ctx.training_dir / "ckpts").is_dir() else []:
        rel = str(path.relative_to(ctx.run_root))
        stale.append({"path": rel, "exists": True, "reason": "downstream_of_ga"})
    ws_root = resolve_workspace_root(workspace_root)
    results_root = ws_root / "model_after_results"
    run_s = str(ctx.run_root.resolve())
    if results_root.is_dir():
        for manifest in results_root.rglob("evaluate_model_manifest.json"):
            try:
                data = _read_json(manifest)
            except (json.JSONDecodeError, OSError):
                continue
            if data.get("model_run_root") == run_s:
                stale.append(
                    {
                        "path": str(manifest.relative_to(ws_root)),
                        "exists": True,
                        "reason": "model_after_results_for_run",
                    }
                )
    return stale


def _write_stale_report(
    ctx: RunContext,
    *,
    trigger: str,
    previous: Optional[Dict[str, Any]],
    new_binding: Dict[str, Any],
    workspace_root: str | Path | None = None,
) -> Optional[Path]:
    stale_files = _collect_stale_files(ctx, workspace_root=workspace_root)
    if not stale_files and previous is None:
        return None
    if previous and previous.get("checksum") == new_binding.get("checksum") and not stale_files:
        return None
    report = {
        "generated_at": _utc_now(),
        "trigger": trigger,
        "previous_binding": previous,
        "new_binding": {
            "ga_set_id": new_binding.get("ga_set_id"),
            "ga_version_id": new_binding.get("ga_version_id"),
            "checksum": new_binding.get("checksum"),
        },
        "stale_files": stale_files,
        "recommended_actions": [
            "python -m app.services.ga_registry_manager apply-ga-to-run --run-root <run>",
            "python -m app.services.orchestrator run-pipeline --run-root <run> --through step13",
        ],
    }
    out = ctx.manifest_dir / "ga_stale_report.json"
    _write_json(out, report)
    return out


def binding_path(ctx: RunContext) -> Path:
    return ctx.manifest_dir / "ga_binding.json"


def materialized_ga_path(ctx: RunContext) -> Path:
    return ctx.run_root / MATERIALIZED_GA_REL


def load_binding(ctx: RunContext) -> Dict[str, Any]:
    path = binding_path(ctx)
    if not path.is_file():
        raise FileNotFoundError(
            f"missing {path.relative_to(ctx.run_root)}; "
            "bind a GA version with: python -m app.services.ga_registry_manager bind-ga-version-to-run"
        )
    return _read_json(path)


def verify_bound_ga(ctx: RunContext) -> Dict[str, Any]:
    binding = load_binding(ctx)
    ga_path = materialized_ga_path(ctx)
    if not ga_path.is_file():
        raise FileNotFoundError(
            f"missing materialized GA at {MATERIALIZED_GA_REL}; "
            "run bind-ga-version-to-run before the core pipeline."
        )
    checksum, num_ga = ga_csv_checksum(ga_path)
    expected = binding.get("checksum")
    if expected and checksum != expected:
        raise ValueError(
            f"GA checksum mismatch: binding has {expected}, "
            f"materialized {ga_path} has {checksum}"
        )
    if binding.get("num_ga") and int(binding["num_ga"]) != num_ga:
        raise ValueError(
            f"GA row count mismatch: binding num_ga={binding['num_ga']}, file has {num_ga}"
        )
    return {
        "ok": True,
        "ga_set_id": binding.get("ga_set_id"),
        "ga_version_id": binding.get("ga_version_id"),
        "checksum": checksum,
        "num_ga": num_ga,
        "materialized_path": str(ga_path.resolve()),
    }


def require_bound_ga(ctx: RunContext) -> Dict[str, Any]:
    try:
        return verify_bound_ga(ctx)
    except (FileNotFoundError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}


def generate_ga_from_run(
    run_root: str | Path,
    *,
    ga_set_id: Optional[str] = None,
    ga_version_id: Optional[str] = None,
    bind_to_run: bool = False,
    workers: Optional[int] = None,
    limit: int = -1,
    workspace_root: str | Path | None = None,
) -> Dict[str, Any]:
    ctx = RunContext.from_run_root(run_root)
    repaired = ctx.tmp_dir / "repaired_ligand_data.csv"
    if not repaired.is_file():
        raise FileNotFoundError(f"missing {repaired}")

    gid = ga_set_id or f"{ctx.run_root.name}_derived"
    if ga_version_id:
        if version_dir(gid, ga_version_id, workspace_root=workspace_root).exists():
            raise FileExistsError(f"version exists: {gid}/{ga_version_id}")
        vid = ga_version_id
    else:
        vid = _next_version_id(gid, workspace_root=workspace_root)

    work = ctx.tmp_dir / "_ga_generate_work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    _run_step4_5(mode="generate-ga", input_dir=ctx.tmp_dir, output_dir=work, workers=workers, limit=limit)
    gen_csv = work / "GA_with_id.csv"
    if not gen_csv.is_file():
        raise RuntimeError("generate-ga did not produce GA_with_id.csv")

    with open(gen_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    rows_norm = _normalize_ga_rows(rows, generate_missing_id=True)
    csv_path = _write_ga_version_files(
        gid,
        vid,
        rows_norm,
        parent_version_id=None,
        source={"type": "from_run_generate_ga", "run_root": str(ctx.run_root.resolve())},
        notes="generate-ga-from-run",
        workspace_root=workspace_root,
    )
    _ensure_ga_set_json(
        gid,
        source={"type": "from_run", "run_root": str(ctx.run_root.resolve())},
        name=gid,
        workspace_root=workspace_root,
    )
    checksum, num_ga = ga_csv_checksum(csv_path)

    result: Dict[str, Any] = {
        "ga_set_id": gid,
        "ga_version_id": vid,
        "ga_csv": str(csv_path.resolve()),
        "checksum": checksum,
        "num_ga": num_ga,
    }
    if bind_to_run:
        result["bind"] = bind_ga_version_to_run(
            ctx.run_root,
            gid,
            vid,
            workspace_root=workspace_root,
        )
    shutil.rmtree(work, ignore_errors=True)
    return result


def create_ga_version_from_csv(
    ga_set_id: str,
    input_csv: str | Path,
    *,
    parent_version_id: Optional[str] = None,
    workspace_root: str | Path | None = None,
) -> Dict[str, Any]:
    sets_root = ga_sets_root(workspace_root)
    if not ga_set_dir(ga_set_id, workspace_root=workspace_root).is_dir() and not (sets_root / ga_set_id).exists():
        _ensure_ga_set_json(ga_set_id, source={"type": "from_csv"}, name=ga_set_id, workspace_root=workspace_root)
    src = Path(input_csv).resolve()
    if not src.is_file():
        raise FileNotFoundError(src)
    with open(src, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "GA_SMILES" not in (reader.fieldnames or []):
            raise ValueError("input CSV must have GA_SMILES column")
        rows = list(reader)
    rows_norm = _normalize_ga_rows(rows, generate_missing_id=True)
    vid = _next_version_id(ga_set_id, workspace_root=workspace_root)
    csv_path = _write_ga_version_files(
        ga_set_id,
        vid,
        rows_norm,
        parent_version_id=parent_version_id,
        source={"type": "from_csv", "input_csv": str(src)},
        notes="create-ga-version-from-csv",
        workspace_root=workspace_root,
    )
    _ensure_ga_set_json(ga_set_id, source={"type": "from_csv"}, name=ga_set_id, workspace_root=workspace_root)
    checksum, num_ga = ga_csv_checksum(csv_path)
    return {
        "ga_set_id": ga_set_id,
        "ga_version_id": vid,
        "ga_csv": str(csv_path.resolve()),
        "checksum": checksum,
        "num_ga": num_ga,
    }


def create_ga_set_or_version_from_upload(
    ga_set_id: str,
    upload_content: bytes | str,
    *,
    parent_version_id: Optional[str] = None,
    workspace_root: str | Path | None = None,
) -> Dict[str, Any]:
    """
    Create a GA version from an uploaded CSV (strict mode).

    Strict mode rules:
    - CSV must contain GA_SMILES and GA_ID columns.
    - Every non-empty row must include both GA_SMILES and GA_ID.
    - Missing GA_ID is not auto-generated in this path.
    """
    if isinstance(upload_content, bytes):
        text = upload_content.decode("utf-8")
    else:
        text = upload_content

    reader = csv.DictReader(StringIO(text))
    fieldnames = reader.fieldnames or []
    if "GA_SMILES" not in fieldnames or "GA_ID" not in fieldnames:
        raise ValueError("uploaded CSV must include GA_SMILES and GA_ID columns")

    rows = list(reader)
    for i, row in enumerate(rows, start=2):
        smi = (row.get("GA_SMILES") or "").strip()
        gid = (row.get("GA_ID") or "").strip()
        if not smi and not gid:
            continue
        if not smi or not gid:
            raise ValueError(f"row {i}: GA_SMILES and GA_ID are both required in strict upload mode")

    rows_norm = _normalize_ga_rows(rows, generate_missing_id=False)
    if not ga_set_dir(ga_set_id, workspace_root=workspace_root).is_dir():
        _ensure_ga_set_json(ga_set_id, source={"type": "upload"}, name=ga_set_id, workspace_root=workspace_root)

    vid = _next_version_id(ga_set_id, workspace_root=workspace_root)
    csv_path = _write_ga_version_files(
        ga_set_id,
        vid,
        rows_norm,
        parent_version_id=parent_version_id,
        source={"type": "upload"},
        notes="create-ga-version-from-upload",
        workspace_root=workspace_root,
    )
    _ensure_ga_set_json(ga_set_id, source={"type": "upload"}, name=ga_set_id, workspace_root=workspace_root)
    checksum, num_ga = ga_csv_checksum(csv_path)
    return {
        "ga_set_id": ga_set_id,
        "ga_version_id": vid,
        "ga_csv": str(csv_path.resolve()),
        "checksum": checksum,
        "num_ga": num_ga,
    }


def bind_ga_version_to_run(
    run_root: str | Path,
    ga_set_id: str,
    ga_version_id: str,
    *,
    bound_by: str = "ga_registry_manager",
    workspace_root: str | Path | None = None,
) -> Dict[str, Any]:
    ctx = RunContext.from_run_root(run_root)
    src_csv = version_dir(ga_set_id, ga_version_id, workspace_root=workspace_root) / "GA_with_id.csv"
    if not src_csv.is_file():
        raise FileNotFoundError(f"GA version not found: {src_csv}")

    previous: Optional[Dict[str, Any]] = None
    bind_file = binding_path(ctx)
    if bind_file.is_file():
        previous = _read_json(bind_file)

    dest = materialized_ga_path(ctx)
    materialize_ga_csv(src_csv, dest)
    checksum, num_ga = ga_csv_checksum(dest)

    payload = {
        "run_root": str(ctx.run_root.resolve()),
        "ga_set_id": ga_set_id,
        "ga_version_id": ga_version_id,
        "source_ga_csv": str(src_csv.resolve()),
        "materialized_path": MATERIALIZED_GA_REL,
        "checksum": checksum,
        "num_ga": num_ga,
        "bound_at": _utc_now(),
        "bound_by": bound_by,
        "status": "bound",
    }
    _write_json(bind_file, payload)

    stale_path = _write_stale_report(
        ctx,
        trigger="bind-ga-version-to-run",
        previous=previous,
        new_binding=payload,
        workspace_root=workspace_root,
    )
    out = dict(payload)
    if stale_path:
        out["stale_report"] = str(stale_path.resolve())
    return out


def apply_ga_to_run(
    run_root: str | Path,
    *,
    workers: Optional[int] = None,
    limit: int = -1,
) -> Dict[str, Any]:
    ctx = RunContext.from_run_root(run_root)
    verify_bound_ga(ctx)
    repaired = ctx.tmp_dir / "repaired_ligand_data.csv"
    if not repaired.is_file():
        raise FileNotFoundError(f"missing {repaired}")
    _run_step4_5(
        mode="apply-ga",
        input_dir=ctx.tmp_dir,
        output_dir=ctx.tmp_dir,
        workers=workers,
        limit=limit,
    )
    for name in ("ligand_with_gac.csv", "IRL_filtered.csv"):
        if not (ctx.tmp_dir / name).is_file():
            raise RuntimeError(f"apply-ga did not produce tmp/{name}")
    return {"ok": True, "run_root": str(ctx.run_root.resolve())}


def _add_workspace_root_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--workspace-root",
        default=None,
        help="Workspace root (default: repo workspace/). GA sets live under {workspace-root}/ga_sets/",
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ChemDB workspace GA registry (Phase 1E)")
    sub = p.add_subparsers(dest="command", required=True)

    g = sub.add_parser("generate-ga-from-run", help="Create workspace GA version from run ligands")
    g.add_argument("--run-root", required=True)
    g.add_argument("--ga-set-id", default=None)
    g.add_argument("--ga-version-id", default=None)
    g.add_argument("--bind-to-run", action="store_true")
    g.add_argument("--workers", type=int, default=None)
    g.add_argument("--limit", type=int, default=-1)
    _add_workspace_root_arg(g)

    c = sub.add_parser("create-ga-version-from-csv", help="New immutable GA version from CSV")
    c.add_argument("--ga-set-id", required=True)
    c.add_argument("--input-csv", required=True)
    c.add_argument("--parent-version-id", default=None)
    _add_workspace_root_arg(c)

    b = sub.add_parser("bind-ga-version-to-run", help="Materialize GA version into run/tmp")
    b.add_argument("--run-root", required=True)
    b.add_argument("--ga-set-id", required=True)
    b.add_argument("--ga-version-id", required=True)
    _add_workspace_root_arg(b)

    a = sub.add_parser("apply-ga-to-run", help="Run step4_5 apply-ga on bound run")
    a.add_argument("--run-root", required=True)
    a.add_argument("--workers", type=int, default=None)
    a.add_argument("--limit", type=int, default=-1)
    _add_workspace_root_arg(a)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    ws = getattr(args, "workspace_root", None)
    try:
        if args.command == "generate-ga-from-run":
            out = generate_ga_from_run(
                args.run_root,
                ga_set_id=args.ga_set_id,
                ga_version_id=args.ga_version_id,
                bind_to_run=args.bind_to_run,
                workers=args.workers,
                limit=args.limit,
                workspace_root=ws,
            )
        elif args.command == "create-ga-version-from-csv":
            out = create_ga_version_from_csv(
                args.ga_set_id,
                args.input_csv,
                parent_version_id=args.parent_version_id,
                workspace_root=ws,
            )
        elif args.command == "bind-ga-version-to-run":
            out = bind_ga_version_to_run(
                args.run_root,
                args.ga_set_id,
                args.ga_version_id,
                workspace_root=ws,
            )
        elif args.command == "apply-ga-to-run":
            out = apply_ga_to_run(args.run_root, workers=args.workers, limit=args.limit)
        else:
            raise SystemExit(f"unknown command: {args.command}")
        print(json.dumps(out, indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
