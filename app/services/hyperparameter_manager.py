"""Hyperparameter resource manager (Phase 3B minimal)."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None  # type: ignore

HP_SET_ID_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{2,63}$")
HP_VERSION_RE = re.compile(r"^v\d{3}$")
WHITELIST_KEYS = {
    "batch_size",
    "epochs",
    "lr",
    "learning_rate",
    "optimizer.lr",
    "optimizer.weight_decay",
    "model.hidden_dim",
    "model.dropout",
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_yaml() -> None:
    if yaml is None:
        raise RuntimeError("pyyaml is required for hyperparameter manager")


def _hp_set_dir(workspace_root: str | Path, hp_set_id: str) -> Path:
    if not HP_SET_ID_RE.match(hp_set_id):
        raise ValueError(f"invalid hp_set_id: {hp_set_id!r}")
    return Path(workspace_root).resolve() / "hyperparameter_sets" / hp_set_id


def _next_version_id(hp_set_dir: Path) -> str:
    versions = hp_set_dir / "versions"
    if not versions.is_dir():
        return "v001"
    nums = []
    for p in versions.iterdir():
        if p.is_dir() and HP_VERSION_RE.match(p.name):
            nums.append(int(p.name[1:]))
    return f"v{max(nums, default=0)+1:03d}"


def create_hyperparameter_version_from_upload(
    *,
    workspace_root: str | Path,
    hp_set_id: str,
    yaml_content: bytes | str,
    notes: str = "",
) -> Dict[str, Any]:
    _require_yaml()
    text = yaml_content.decode("utf-8") if isinstance(yaml_content, bytes) else yaml_content
    data = yaml.safe_load(text)  # type: ignore[attr-defined]
    if not isinstance(data, dict):
        raise ValueError("hyperparameter yaml must be a mapping")
    sdir = _hp_set_dir(workspace_root, hp_set_id)
    sdir.mkdir(parents=True, exist_ok=True)
    meta = sdir / "hp_set.json"
    now = _utc_now()
    if not meta.is_file():
        meta.write_text(
            json.dumps(
                {
                    "hp_set_id": hp_set_id,
                    "name": hp_set_id,
                    "description": "",
                    "created_at": now,
                    "updated_at": now,
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
    else:
        payload = json.loads(meta.read_text(encoding="utf-8"))
        payload["updated_at"] = now
        meta.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    vid = _next_version_id(sdir)
    vdir = sdir / "versions" / vid
    vdir.mkdir(parents=True, exist_ok=False)
    ypath = vdir / "config.yaml"
    ypath.write_text(text, encoding="utf-8")
    (vdir / "hp_version.json").write_text(
        json.dumps(
            {
                "hp_set_id": hp_set_id,
                "hp_version_id": vid,
                "created_at": now,
                "created_by": "hyperparameter_manager",
                "notes": notes,
                "source": {"type": "upload"},
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return {
        "hp_set_id": hp_set_id,
        "hp_version_id": vid,
        "yaml_path": str(ypath.resolve()),
    }


def _flatten(prefix: str, obj: Any) -> List[Tuple[str, Any]]:
    if isinstance(obj, dict):
        out: List[Tuple[str, Any]] = []
        for k, v in obj.items():
            path = f"{prefix}.{k}" if prefix else str(k)
            out.extend(_flatten(path, v))
        return out
    return [(prefix, obj)]


def _set_nested(cfg: Dict[str, Any], path: str, value: Any) -> bool:
    parts = path.split(".")
    cur: Any = cfg
    for p in parts[:-1]:
        if not isinstance(cur, dict) or p not in cur:
            return False
        cur = cur[p]
    if not isinstance(cur, dict):
        return False
    cur[parts[-1]] = value
    return True


def apply_hyperparameters_to_run(
    *,
    workspace_root: str | Path,
    run_root: str | Path,
    hp_set_id: str,
    hp_version_id: str,
) -> Dict[str, Any]:
    _require_yaml()
    ws = Path(workspace_root).resolve()
    rroot = Path(run_root).resolve()
    hp_yaml = ws / "hyperparameter_sets" / hp_set_id / "versions" / hp_version_id / "config.yaml"
    if not hp_yaml.is_file():
        raise FileNotFoundError(f"hyperparameter yaml not found: {hp_yaml}")
    cfg_path = rroot / "training" / "config.yaml"
    if not cfg_path.is_file():
        raise FileNotFoundError(f"training config not found: {cfg_path}")

    hp = yaml.safe_load(hp_yaml.read_text(encoding="utf-8"))  # type: ignore[attr-defined]
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))  # type: ignore[attr-defined]
    if not isinstance(hp, dict) or not isinstance(cfg, dict):
        raise ValueError("invalid yaml structure")

    applied: List[str] = []
    ignored: List[str] = []
    for key, value in _flatten("", hp):
        if key in WHITELIST_KEYS and _set_nested(cfg, key, value):
            applied.append(key)
        else:
            ignored.append(key)

    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")  # type: ignore[attr-defined]
    manifest = {
        "hp_set_id": hp_set_id,
        "hp_version_id": hp_version_id,
        "applied": applied,
        "ignored": ignored,
        "applied_at": _utc_now(),
        "source_yaml": str(hp_yaml.resolve()),
    }
    out = rroot / "manifests" / "hyperparameter_binding.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Hyperparameter manager")
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("create-hp-version-from-yaml")
    c.add_argument("--workspace-root", required=True)
    c.add_argument("--hp-set-id", required=True)
    c.add_argument("--input-yaml", required=True)
    c.add_argument("--notes", default="")
    a = sub.add_parser("apply-to-run")
    a.add_argument("--workspace-root", required=True)
    a.add_argument("--run-root", required=True)
    a.add_argument("--hp-set-id", required=True)
    a.add_argument("--hp-version-id", required=True)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "create-hp-version-from-yaml":
        content = Path(args.input_yaml).read_bytes()
        out = create_hyperparameter_version_from_upload(
            workspace_root=args.workspace_root,
            hp_set_id=args.hp_set_id,
            yaml_content=content,
            notes=args.notes,
        )
        print(json.dumps(out, ensure_ascii=False, indent=2))
        return 0
    if args.command == "apply-to-run":
        out = apply_hyperparameters_to_run(
            workspace_root=args.workspace_root,
            run_root=args.run_root,
            hp_set_id=args.hp_set_id,
            hp_version_id=args.hp_version_id,
        )
        print(json.dumps(out, ensure_ascii=False, indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
