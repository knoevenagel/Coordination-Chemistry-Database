from __future__ import annotations

import streamlit as st

from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Model Library - ChemDB Web", layout="wide")
render_sidebar()
st.title("Model Library")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found. Go to Settings and rebuild index.")
else:
    try:
        all_runs = repositories.list_runs(conn)
        project_options = sorted({r["project_id"] for r in all_runs})
    finally:
        conn.close()

    project_id = st.selectbox("Filter by project_id", ["(all)"] + project_options)
    run_options = sorted({r["run_id"] for r in all_runs if project_id == "(all)" or r["project_id"] == project_id})
    run_id = st.selectbox("Filter by run_id", ["(all)"] + run_options)

    q_project = None if project_id == "(all)" else project_id
    q_run = None if run_id == "(all)" else run_id

    conn = state.get_db_connection(read_only=True)
    if conn is None:
        st.warning("SQLite DB disconnected. Refresh page.")
    else:
        try:
            models = repositories.list_models(conn, project_id=q_project, run_id=q_run)
        finally:
            conn.close()

        st.subheader("Models")
        show_table(
            [
                {
                    "model_id": r.get("model_id"),
                    "project_id": r.get("project_id"),
                    "run_id": r.get("run_id"),
                    "checkpoint_stem": r.get("checkpoint_stem"),
                    "embedding_backend": r.get("embedding_backend"),
                    "val_mrr": r.get("val_mrr"),
                    "val_recall_at_1": r.get("val_recall_at_1"),
                    "val_recall_at_5": r.get("val_recall_at_5"),
                    "ga_set_id": r.get("ga_set_id"),
                    "ga_version_id": r.get("ga_version_id"),
                    "size_bytes": r.get("size_bytes"),
                    "ingest_status": r.get("ingest_status"),
                }
                for r in models
            ],
            empty_message="No models found.",
        )

        if models:
            model_id = st.selectbox("Select Model", [m["model_id"] for m in models])
            conn = state.get_db_connection(read_only=True)
            if conn is None:
                st.warning("SQLite DB disconnected. Refresh page.")
            else:
                try:
                    model = repositories.get_model(conn, model_id)
                finally:
                    conn.close()

                st.subheader("Model Detail")
                st.json(
                    {
                        "model_id": model.get("model_id") if model else None,
                        "checkpoint_path": model.get("checkpoint_path") if model else None,
                        "config_path": model.get("config_path") if model else None,
                        "index_path": model.get("index_path") if model else None,
                        "l3_embedding_path": model.get("l3_embedding_path") if model else None,
                        "metal_embedding_path": model.get("metal_embedding_path") if model else None,
                        "history_path": model.get("history_path") if model else None,
                        "ga_binding_checksum": model.get("ga_binding_checksum") if model else None,
                        "ingest_error": model.get("ingest_error") if model else None,
                    }
                )
