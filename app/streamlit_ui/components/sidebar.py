from __future__ import annotations

import sqlite3

import streamlit as st

from app.storage import repositories
from app.streamlit_ui import state


def render_sidebar() -> None:
    status = state.get_db_status()
    with st.sidebar:
        st.title("ChemDB Web")
        st.caption("Phase 3B MVP")
        st.markdown(f"**Workspace:** `{status['workspace_root']}`")
        st.markdown(f"**DB:** `{status['db_path']}`")
        if status["db_exists"]:
            st.success("SQLite DB exists")
            conn = state.get_db_connection(read_only=True)
            if conn is not None:
                try:
                    try:
                        recent_jobs = repositories.get_recent_jobs(conn, limit=3)
                    except sqlite3.OperationalError:
                        recent_jobs = []
                finally:
                    conn.close()
                if recent_jobs:
                    st.markdown("**Recent Jobs:**")
                    for job in recent_jobs:
                        st.caption(f"{job.get('status')} · {job.get('title') or job.get('job_type')} ({job.get('job_id')})")
        else:
            st.warning("SQLite DB missing. Go to Settings to rebuild index.")
