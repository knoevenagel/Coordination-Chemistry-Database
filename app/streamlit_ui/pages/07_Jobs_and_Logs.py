from __future__ import annotations

import sqlite3
from pathlib import Path

import streamlit as st

from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Jobs & Logs - ChemDB Web", layout="wide")
render_sidebar()
st.title("Jobs & Logs")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found. Go to Settings and rebuild index.")
else:
    try:
        try:
            jobs = repositories.get_recent_jobs(conn, limit=50)
        except sqlite3.OperationalError:
            jobs = []
    finally:
        conn.close()

    show_table(
        [
            {
                "job_id": j.get("job_id"),
                "job_type": j.get("job_type"),
                "title": j.get("title"),
                "status": j.get("status"),
                "created_at": j.get("created_at"),
                "started_at": j.get("started_at"),
                "finished_at": j.get("finished_at"),
                "log_path": j.get("log_path"),
                "error": j.get("error"),
            }
            for j in jobs
        ],
        empty_message="No jobs found.",
    )

    if jobs:
        labels = [f"{j['job_id']} · {j.get('status')} · {j.get('title') or j.get('job_type')}" for j in jobs]
        idx = st.selectbox("Select Job", list(range(len(labels))), format_func=lambda i: labels[i])
        job = jobs[idx]
        st.subheader("Job Detail")
        st.json(job)

        st.subheader("Log Preview")
        log_path = job.get("log_path")
        if not log_path:
            st.warning("No log path recorded.")
        else:
            p = Path(log_path)
            if not p.is_file():
                st.warning(f"log file not found: {p}")
            else:
                try:
                    text = p.read_text(encoding="utf-8")
                except Exception as exc:  # pragma: no cover
                    st.warning(str(exc))
                else:
                    lines = text.splitlines()
                    tail_n = st.slider("Lines from end", min_value=20, max_value=1000, value=200, step=20)
                    st.code("\n".join(lines[-tail_n:]), language="text")
