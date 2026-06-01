from __future__ import annotations

import subprocess
import sys

import streamlit as st

from app.services import task_registry_manager
from app.storage import previews, repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.previews import show_csv_preview
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Task Library - ChemDB Web", layout="wide")
render_sidebar()
st.title("Task Library")

status = state.get_db_status()
with st.expander("Create Task From Upload", expanded=False):
    with st.form("task_upload_form"):
        new_task_id = st.text_input("Task ID", value="")
        metal = st.text_input("Metal", value="")
        embedding = st.text_input("Embedding", value="")
        notes = st.text_area("Notes", value="")
        candidates_file = st.file_uploader("Candidates CSV", type=["csv"])
        positives_file = st.file_uploader("Positives CSV", type=["csv"])
        negatives_file = st.file_uploader("Negatives CSV (optional)", type=["csv"])
        submit_task = st.form_submit_button("Create Task")
    if submit_task:
        if not new_task_id.strip():
            st.error("Task ID is required.")
        elif candidates_file is None or positives_file is None:
            st.error("Candidates and positives CSV are required.")
        else:
            try:
                result = task_registry_manager.create_task_from_uploads(
                    task_id=new_task_id.strip(),
                    workspace_root=status["workspace_root"],
                    candidates_csv=candidates_file.getvalue(),
                    positives_csv=positives_file.getvalue(),
                    negatives_csv=negatives_file.getvalue() if negatives_file is not None else None,
                    metal=metal.strip(),
                    embedding=embedding.strip(),
                    notes=notes.strip(),
                )
                st.success(
                    f"Created task {result['task_id']} "
                    f"(candidates={result['candidate_count']}, positives={result['positive_count']}, "
                    f"negatives={result['negative_count']})"
                )
            except Exception as exc:
                st.error(f"Create task failed: {exc}")
            else:
                cmd = [
                    sys.executable,
                    "-m",
                    "app.storage.cli",
                    "rebuild-index",
                    "--workspace-root",
                    status["workspace_root"],
                    "--db-path",
                    status["db_path"],
                ]
                with st.spinner("Rebuilding SQLite index..."):
                    proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode == 0:
                    st.success("Rebuild index completed.")
                else:
                    st.error(f"Rebuild index failed (exit={proc.returncode})")
                    st.code(proc.stderr or "(empty)", language="text")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found. Go to Settings and rebuild index.")
else:
    try:
        tasks = repositories.list_model_after_tasks(conn)
    finally:
        conn.close()

    st.subheader("Tasks")
    show_table(
        [
            {
                "task_id": r.get("task_id"),
                "metal": r.get("metal"),
                "embedding": r.get("embedding"),
                "candidate_count": r.get("candidate_count"),
                "positive_count": r.get("positive_count"),
                "negative_count": r.get("negative_count"),
                "notes": r.get("notes"),
                "ingest_status": r.get("ingest_status"),
            }
            for r in tasks
        ],
        empty_message="No tasks found.",
    )

    if tasks:
        task_id = st.selectbox("Select Task", [r["task_id"] for r in tasks])
        conn = state.get_db_connection(read_only=True)
        if conn is None:
            st.warning("SQLite DB disconnected. Refresh page.")
        else:
            try:
                task = repositories.get_model_after_task(conn, task_id)
                batches = repositories.list_batches(conn, task_id)
            finally:
                conn.close()

            st.subheader("Task Detail")
            st.json(task or {})

            if task:
                task_previews = previews.preview_task_files(task, limit=100)
                tab_candidates, tab_positives, tab_negatives = st.tabs(["Candidates", "Positives", "Negatives"])
                with tab_candidates:
                    show_csv_preview(task_previews["candidates"], "Candidates Preview")
                with tab_positives:
                    show_csv_preview(task_previews["positives"], "Positives Preview")
                with tab_negatives:
                    show_csv_preview(task_previews["negatives"], "Negatives Preview")

            st.subheader("Task Batches")
            show_table(batches, empty_message="No batches for this task.")
