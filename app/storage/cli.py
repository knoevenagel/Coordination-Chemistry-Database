"""Phase 2A CLI: read-only workspace → SQLite registry."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from .db import init_db
from . import ingest
from . import repositories

WEB_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORKSPACE = WEB_ROOT / "workspace"
DEFAULT_DB = DEFAULT_WORKSPACE / "chemdb.sqlite"


def _print_json(payload: Dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def cmd_init_db(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    init_db(db_path)
    _print_json({"ok": True, "db_path": str(db_path)})
    return 0


def cmd_register_project(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    ws = Path(args.workspace_root).resolve()
    conn = init_db(db_path)
    result = ingest.ingest_project(conn, ws, args.project_id)
    conn.close()
    _print_json(result)
    return 0 if result.get("ingest_status") != "error" else 1


def cmd_register_run(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    run_root = Path(args.run_root).resolve()
    conn = init_db(db_path)
    result = ingest.register_run(conn, run_root)
    conn.close()
    _print_json(result)
    return 0


def cmd_register_ga_set(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    ws = Path(args.workspace_root).resolve()
    conn = init_db(db_path)
    result = ingest.ingest_ga_set(conn, ws, args.ga_set_id)
    conn.close()
    _print_json(result)
    return 0 if result.get("ingest_status") != "error" else 1


def cmd_register_task(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    task_dir = Path(args.task_dir).resolve()
    conn = init_db(db_path)
    result = ingest.ingest_model_after_task(conn, task_dir)
    conn.close()
    _print_json(result)
    return 0


def cmd_register_model_after_batch(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    batch_dir = Path(args.batch_dir).resolve()
    conn = init_db(db_path)
    result = ingest.ingest_model_after_batch(conn, batch_dir)
    conn.close()
    _print_json(result)
    return 0


def cmd_rebuild_index(args: argparse.Namespace) -> int:
    db_path = Path(args.db_path).resolve()
    ws = Path(args.workspace_root).resolve()
    conn = init_db(db_path)
    summary = ingest.rebuild_index(conn, ws)
    conn.close()
    _print_json({"ok": True, "summary": summary})
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ChemDB Phase 2A SQLite registry (read-only ingest)")
    sub = p.add_subparsers(dest="command", required=True)

    db_kw = {"default": str(DEFAULT_DB), "help": "SQLite database path"}
    ws_kw = {"default": str(DEFAULT_WORKSPACE), "help": "Workspace root directory"}

    s = sub.add_parser("init-db", help="Create schema in chemdb.sqlite")
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_init_db)

    s = sub.add_parser("register-project", help="Ingest one project")
    s.add_argument("--workspace-root", **ws_kw)
    s.add_argument("--project-id", required=True)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_register_project)

    s = sub.add_parser("register-run", help="Ingest one run (full run bundle)")
    s.add_argument("--run-root", required=True)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_register_run)

    s = sub.add_parser("register-ga-set", help="Ingest GA set and all versions")
    s.add_argument("--workspace-root", **ws_kw)
    s.add_argument("--ga-set-id", required=True)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_register_ga_set)

    s = sub.add_parser("register-task", help="Ingest model-after task")
    s.add_argument("--task-dir", required=True)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_register_task)

    s = sub.add_parser("register-model-after-batch", help="Ingest model-after batch results")
    s.add_argument("--batch-dir", required=True)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_register_model_after_batch)

    s = sub.add_parser("rebuild-index", help="Full workspace scan and UPSERT")
    s.add_argument("--workspace-root", **ws_kw)
    s.add_argument("--db-path", **db_kw)
    s.set_defaults(func=cmd_rebuild_index)

    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
