from __future__ import annotations

import streamlit as st

from app.storage import previews, repositories
from app.streamlit_ui import state
from app.streamlit_ui.components.previews import show_csv_preview, show_json_preview
from app.streamlit_ui.components.sidebar import render_sidebar
from app.streamlit_ui.components.tables import show_table

st.set_page_config(page_title="Model Selector Results - ChemDB Web", layout="wide")
render_sidebar()
st.title("Model Selector Results")

conn = state.get_db_connection(read_only=True)
if conn is None:
    st.warning("SQLite DB not found. Go to Settings and rebuild index.")
else:
    try:
        tasks = repositories.list_model_after_tasks(conn)
    finally:
        conn.close()

    if not tasks:
        st.info("No model-after tasks found.")
    else:
        task_id = st.selectbox("Select Task", [t["task_id"] for t in tasks])
        conn = state.get_db_connection(read_only=True)
        if conn is None:
            st.warning("SQLite DB disconnected. Refresh page.")
        else:
            try:
                batches = repositories.list_batches(conn, task_id)
            finally:
                conn.close()

            if not batches:
                st.info("No batches for this task.")
            else:
                batch_id = st.selectbox("Select Batch", [b["batch_id"] for b in batches])
                conn = state.get_db_connection(read_only=True)
                if conn is None:
                    st.warning("SQLite DB disconnected. Refresh page.")
                else:
                    try:
                        batch = repositories.get_batch(conn, task_id, batch_id)
                        results = repositories.list_model_results_with_models(conn, task_id, batch_id)
                    finally:
                        conn.close()

                    st.subheader("Batch Metadata")
                    st.json(
                        {
                            "command": batch.get("command") if batch else None,
                            "finished_at": batch.get("finished_at") if batch else None,
                            "n_models": batch.get("n_models") if batch else None,
                            "n_success": batch.get("n_success") if batch else None,
                            "best_model_id": batch.get("best_model_id") if batch else None,
                            "output_dir": batch.get("output_dir") if batch else None,
                            "summary_path": batch.get("summary_path") if batch else None,
                        }
                    )

                    st.subheader("Model Ranking")
                    show_table(
                        [
                            {
                                "rank_among_models": r.get("rank_among_models"),
                                "raw_model_id": r.get("raw_model_id"),
                                "registry_model_id": r.get("registry_model_id"),
                                "model_join_status": r.get("model_join_status"),
                                "model_name": r.get("model_name"),
                                "project_id": r.get("project_id"),
                                "run_id": r.get("run_id"),
                                "checkpoint_stem": r.get("checkpoint_stem"),
                                "ga_set_id": r.get("ga_set_id"),
                                "ga_version_id": r.get("ga_version_id"),
                                "mrr": r.get("mrr"),
                                "hit_at_5": r.get("hit_at_5"),
                                "hit_at_10": r.get("hit_at_10"),
                                "hit_at_20": r.get("hit_at_20"),
                                "hit_at_50": r.get("hit_at_50"),
                                "selection_score": r.get("selection_score"),
                                "status": r.get("status"),
                            }
                            for r in results
                        ],
                        empty_message="No model results found.",
                    )

                    if results:
                        choice_label = [f"{r.get('raw_model_id')} (rank={r.get('rank_among_models')})" for r in results]
                        choice_idx = st.selectbox(
                            "Select Result",
                            list(range(len(choice_label))),
                            format_func=lambda i: choice_label[i],
                        )
                        selected = results[choice_idx]
                        pv = previews.preview_model_result(selected, limit=100)

                        tab_ranking, tab_metrics, tab_report = st.tabs(["ranking.csv", "metrics.json", "report.json"])
                        with tab_ranking:
                            show_csv_preview(pv["ranking"], "Ranking Preview (first 100 rows)")
                        with tab_metrics:
                            show_json_preview(pv["metrics"], "Metrics JSON")
                        with tab_report:
                            show_json_preview(pv["report"], "Report JSON")
