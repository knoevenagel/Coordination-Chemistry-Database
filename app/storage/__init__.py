"""Phase 2A SQLite read-only workspace registry."""

from .db import connect_db, execute_schema, init_db

__all__ = ["connect_db", "execute_schema", "init_db"]
