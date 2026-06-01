"""Shared UI state and DB access helpers."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from typing import Optional

from app.storage.db import connect_db


def get_workspace_root() -> Path:
    raw = os.getenv("CHEMDB_WORKSPACE_ROOT", "workspace")
    return Path(raw).expanduser().resolve()


def get_db_path() -> Path:
    raw = os.getenv("CHEMDB_DB_PATH")
    if raw:
        return Path(raw).expanduser().resolve()
    return (get_workspace_root() / "chemdb.sqlite").resolve()


def get_db_status() -> dict:
    workspace_root = get_workspace_root()
    db_path = get_db_path()
    return {
        "workspace_root": str(workspace_root),
        "db_path": str(db_path),
        "db_exists": db_path.is_file(),
    }


def get_db_connection(read_only: bool = True) -> Optional[sqlite3.Connection]:
    db_path = get_db_path()
    if read_only and not db_path.is_file():
        return None
    try:
        return connect_db(db_path, read_only=read_only)
    except sqlite3.Error:
        return None
