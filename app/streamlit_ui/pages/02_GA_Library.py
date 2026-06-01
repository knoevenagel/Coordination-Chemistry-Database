from __future__ import annotations

import subprocess
import sys

import streamlit as st

from app.services import ga_registry_manager
from app.storage import previews, repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.previews import show_csv_preview
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="GA Library - ChemDB Web", layout="wide")
render_sidebar()
st.title("GA Library")

status = state.get_db_status()
with st.expander("Create GA Version From Upload (strict GA_SMILES + GA_ID)", expanded=False):
    with st.form("ga_upload_form"):
        new_ga_set_id = st.text_input("GA Set ID", value="")
        uploaded_csv = st.file_uploader("Upload GA CSV", type=["csv"])
        submit_upload = st.form_submit_button("Create GA Version")
    if submit_upload:
        if not new_ga_set_id.strip():
            st.error("GA Set ID is required.")
        elif uploaded_csv is None:
            st.error("Please upload a CSV file.")
        else:
            try:
                result = ga_registry_manager.create_ga_set_or_version_from_upload(
                    new_ga_set_id.strip(),
                    uploaded_csv.getvalue(),
                    workspace_root=status["workspace_root"],
                )
                st.success(
                    f"Created {result['ga_set_id']}/{result['ga_version_id']} "
                    f"(num_ga={result['num_ga']}, checksum={result['checksum']})"
                )
            except Exception as exc:
                st.error(f"Create GA version failed: {exc}")
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
        ga_sets = repositories.list_ga_sets(conn)
    finally:
        conn.close()

    st.subheader("GA Sets")
    show_table(
        [
            {
                "ga_set_id": r.get("ga_set_id"),
                "name": r.get("name"),
                "description": r.get("description"),
                "created_at": r.get("created_at"),
                "updated_at": r.get("updated_at"),
                "ingest_status": r.get("ingest_status"),
            }
            for r in ga_sets
        ],
        empty_message="No GA sets found.",
    )

    if ga_sets:
        ga_set_id = st.selectbox("Select GA Set", [r["ga_set_id"] for r in ga_sets])
        conn = state.get_db_connection(read_only=True)
        if conn is None:
            st.warning("SQLite DB disconnected. Refresh page.")
        else:
            try:
                versions = repositories.list_ga_versions(conn, ga_set_id)
            finally:
                conn.close()

            st.subheader("GA Versions")
            show_table(
                [
                    {
                        "ga_version_id": r.get("ga_version_id"),
                        "num_ga": r.get("num_ga"),
                        "checksum": r.get("checksum"),
                        "created_at": r.get("created_at"),
                        "created_by": r.get("created_by"),
                        "notes": r.get("notes"),
                        "ingest_status": r.get("ingest_status"),
                    }
                    for r in versions
                ],
                empty_message="No versions found.",
            )

            if versions:
                ga_version_id = st.selectbox("Select GA Version", [r["ga_version_id"] for r in versions])
                conn = state.get_db_connection(read_only=True)
                if conn is None:
                    st.warning("SQLite DB disconnected. Refresh page.")
                else:
                    try:
                        ga_preview = previews.preview_ga_version(conn, ga_set_id, ga_version_id, limit=100)
                        runs_using_version = repositories.list_runs_by_ga_version(conn, ga_set_id, ga_version_id)
                    finally:
                        conn.close()

                    show_csv_preview(ga_preview, "GA CSV Preview (first 100 rows)")
                    st.subheader("Runs Using This GA Version")
                    show_table(runs_using_version, empty_message="No runs bound to this GA version.")
