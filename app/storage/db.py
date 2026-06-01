"""SQLite connection and schema bootstrap for Phase 2A registry."""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]

_SCHEMA_PATH = Path(__file__).resolve().parent / "schema.sql"


def connect_db(db_path: PathLike, *, read_only: bool = False) -> sqlite3.Connection:
    path = Path(db_path).resolve()
    if read_only:
        uri = f"file:{path}?mode=ro"
        conn = sqlite3.connect(uri, uri=True)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def execute_schema(conn: sqlite3.Connection, schema_path: Optional[PathLike] = None) -> None:
    path = Path(schema_path) if schema_path else _SCHEMA_PATH
    sql = path.read_text(encoding="utf-8")
    conn.executescript(sql)
    conn.commit()


def init_db(db_path: PathLike, schema_path: Optional[PathLike] = None) -> sqlite3.Connection:
    conn = connect_db(db_path)
    execute_schema(conn, schema_path)
    return conn


def table_names(conn: sqlite3.Connection) -> list[str]:
    rows = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
    ).fetchall()
    return [r["name"] for r in rows]
