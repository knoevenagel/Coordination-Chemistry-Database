from __future__ import annotations

import streamlit as st

from app.services import model_after_manager
from app.storage import repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Model Selector - ChemDB Web", layout="wide")
render_sidebar()
st.title("Model Selector")
st.caption("Launch evaluate-models job and ingest results back to SQLite.")

status = state.get_db_status()
conn = state.get_db_connection(read_only=True)
tasks = []
models = []
if conn is not None:
    try:
        tasks = repositories.list_model_after_tasks(conn)
        models = repositories.list_models(conn)
    finally:
        conn.close()

if not tasks:
    st.info("No tasks found in registry.")
else:
    task_id = st.selectbox("Task", [t["task_id"] for t in tasks])
    st.subheader("Available Models")
    show_table(
        [
            {
                "model_id": m.get("model_id"),
                "project_id": m.get("project_id"),
                "run_id": m.get("run_id"),
                "checkpoint_stem": m.get("checkpoint_stem"),
                "ga_set_id": m.get("ga_set_id"),
                "ga_version_id": m.get("ga_version_id"),
            }
            for m in models
        ],
        empty_message="No models found.",
    )
    model_options = [m["model_id"] for m in models]
    selected_ids = st.multiselect("Select Models", model_options, default=model_options[: min(5, len(model_options))])
    batch_id = st.text_input("Batch ID (optional)", value="")
    device = st.selectbox("Device", ["cpu", "auto", "cuda"], index=0)

    if st.button("Submit Model Selector Job"):
        selected_models = [m for m in models if m["model_id"] in selected_ids]
        if not selected_models:
            st.error("Select at least one model.")
        else:
            try:
                out = model_after_manager.submit_evaluate_models_job(
                    db_path=status["db_path"],
                    workspace_root=status["workspace_root"],
                    task_id=task_id,
                    models=selected_models,
                    batch_id=(batch_id.strip() or None),
                    device=device,
                )
            except Exception as exc:
                st.error(f"Submit failed: {exc}")
            else:
                st.success(f"Submitted job: {out['job']['job_id']}")
                st.json(out)
                st.info("Track in `07_Jobs_and_Logs`; results appear in `05_Model_Selector_Results` after completion.")
