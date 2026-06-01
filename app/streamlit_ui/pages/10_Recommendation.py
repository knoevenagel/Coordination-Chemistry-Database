from __future__ import annotations

from pathlib import Path

import streamlit as st

from app.services import model_after_manager, task_registry_manager
from app.storage import previews, repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.previews import show_csv_preview
from app.streamlit_ui.components.sidebar import render_sidebar

st.set_page_config(page_title="Recommendation - ChemDB Web", layout="wide")
render_sidebar()
st.title("Recommendation")
st.caption("Run recommend job for existing or temporary task.")

status = state.get_db_status()
ws = Path(status["workspace_root"]).resolve()
conn = state.get_db_connection(read_only=True)
tasks = []
models = []
if conn is not None:
    try:
        tasks = repositories.list_model_after_tasks(conn)
        models = repositories.list_models(conn)
    finally:
        conn.close()

task_mode = st.radio("Task Source", ["Existing Task", "Temporary Upload"], horizontal=True)
task_dir: Path | None = None
cleanup_temp = False
if task_mode == "Existing Task":
    if tasks:
        task_id = st.selectbox("Task", [t["task_id"] for t in tasks])
        task_dir = ws / "model_after_tasks" / task_id
    else:
        st.warning("No tasks available.")
else:
    temp_task_id = st.text_input("Temporary Task ID", value="tmp_task_recommend")
    c_file = st.file_uploader("Candidates CSV", type=["csv"], key="rec_c")
    p_file = st.file_uploader("Positives CSV", type=["csv"], key="rec_p")
    n_file = st.file_uploader("Negatives CSV (optional)", type=["csv"], key="rec_n")
    if st.button("Create Temporary Task"):
        if not temp_task_id.strip() or c_file is None or p_file is None:
            st.error("Temporary task id / candidates / positives are required.")
        else:
            try:
                out = task_registry_manager.create_task_from_uploads(
                    task_id=temp_task_id.strip(),
                    workspace_root=ws,
                    candidates_csv=c_file.getvalue(),
                    positives_csv=p_file.getvalue(),
                    negatives_csv=n_file.getvalue() if n_file else None,
                    notes="temporary recommendation task",
                )
            except FileExistsError:
                task_dir = ws / "model_after_tasks" / temp_task_id.strip()
                st.info("Task already exists, reuse it.")
            except Exception as exc:
                st.error(f"Create temporary task failed: {exc}")
            else:
                task_dir = Path(out["task_dir"])
                cleanup_temp = True
                st.success(f"Temporary task ready: {task_dir}")

if models:
    model_id = st.selectbox("Model", [m["model_id"] for m in models])
    chosen = next((m for m in models if m["model_id"] == model_id), None)
    model_run_root = ws / "projects" / str(chosen.get("project_id")) / "runs" / str(chosen.get("run_id")) if chosen else None
else:
    model_id = ""
    model_run_root = None
    st.warning("No models available in registry.")

batch_id = st.text_input("Batch ID (optional)", value="")
if st.button("Submit Recommendation Job"):
    if task_dir is None or not task_dir.is_dir():
        st.error("Task directory is not ready.")
    elif not model_run_root or not model_run_root.is_dir():
        st.error("Model run root not found.")
    else:
        try:
            out = model_after_manager.submit_recommend_job(
                db_path=status["db_path"],
                workspace_root=ws,
                task_dir=task_dir,
                model_run_root=model_run_root,
                model_id=model_id or None,
                batch_id=(batch_id.strip() or None),
            )
        except Exception as exc:
            st.error(f"Submit failed: {exc}")
        else:
            st.success(f"Submitted job: {out['job']['job_id']}")
            st.json(out)
            if cleanup_temp:
                st.caption("Temporary task can be manually removed after result generation.")

st.subheader("Existing Recommendation Results")
results_root = ws / "recommendation_results"
result_dirs = sorted([p for p in results_root.rglob("*") if p.is_dir() and (p / "ranking.csv").is_file()])
if not result_dirs:
    st.info("No recommendation results found.")
else:
    labels = [str(p.relative_to(results_root)) for p in result_dirs]
    idx = st.selectbox("Select Result Folder", list(range(len(labels))), format_func=lambda i: labels[i])
    result_dir = result_dirs[idx]
    ranking_path = result_dir / "ranking.csv"
    pv = previews.read_csv_preview(ranking_path, limit=100)
    show_csv_preview(pv, "Ranking Preview (first 100 rows)")
    st.download_button(
        "Download ranking.csv",
        data=ranking_path.read_bytes(),
        file_name=f"{result_dir.name}_ranking.csv",
        mime="text/csv",
    )
