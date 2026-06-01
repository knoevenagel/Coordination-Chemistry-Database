-- ChemDBWebVersion Phase 2A: read-only workspace registry (paths + summaries only).

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    project_id TEXT PRIMARY KEY,
    project_root TEXT NOT NULL,
    name TEXT,
    created_at TEXT,
    updated_at TEXT,
    notes TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS runs (
    project_id TEXT NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
    run_id TEXT NOT NULL,
    run_root TEXT NOT NULL,
    created_at TEXT,
    chemdb_repo_root TEXT,
    status TEXT,
    pubchem_file_count INTEGER,
    pipeline_core_ok INTEGER,
    pipeline_training_ok INTEGER,
    last_finished_at TEXT,
    ga_binding_status TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (project_id, run_id)
);

CREATE TABLE IF NOT EXISTS ga_sets (
    ga_set_id TEXT PRIMARY KEY,
    ga_set_root TEXT NOT NULL,
    name TEXT,
    description TEXT,
    created_at TEXT,
    updated_at TEXT,
    tags_json TEXT,
    source_json TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS ga_versions (
    ga_set_id TEXT NOT NULL REFERENCES ga_sets(ga_set_id) ON DELETE CASCADE,
    ga_version_id TEXT NOT NULL,
    parent_version_id TEXT,
    created_at TEXT,
    created_by TEXT,
    num_ga INTEGER,
    checksum TEXT,
    ga_csv_path TEXT,
    source_json TEXT,
    notes TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (ga_set_id, ga_version_id)
);

CREATE TABLE IF NOT EXISTS run_ga_bindings (
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    ga_set_id TEXT,
    ga_version_id TEXT,
    checksum TEXT,
    num_ga INTEGER,
    source_ga_csv TEXT,
    materialized_path TEXT,
    bound_at TEXT,
    bound_by TEXT,
    status TEXT NOT NULL,
    is_current INTEGER NOT NULL DEFAULT 1,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (project_id, run_id),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE,
    FOREIGN KEY (ga_set_id) REFERENCES ga_sets(ga_set_id),
    CHECK (is_current IN (0, 1))
);

CREATE TABLE IF NOT EXISTS step_executions (
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    step_id TEXT NOT NULL,
    started_at TEXT NOT NULL DEFAULT '',
    status TEXT,
    command_json TEXT,
    cwd TEXT,
    finished_at TEXT,
    duration_seconds REAL,
    exit_code INTEGER,
    error TEXT,
    log_path TEXT,
    manifest_path TEXT,
    in_process_result_json TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (project_id, run_id, step_id, started_at),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS artifacts (
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    step_id TEXT NOT NULL DEFAULT '',
    rel_path TEXT NOT NULL,
    role TEXT NOT NULL,
    abs_path TEXT,
    size_bytes INTEGER,
    exists_flag INTEGER NOT NULL DEFAULT 0,
    artifact_type TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (project_id, run_id, step_id, rel_path, role),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE,
    CHECK (role IN ('input', 'output', 'key'))
);

CREATE TABLE IF NOT EXISTS models (
    model_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    checkpoint_path TEXT NOT NULL,
    checkpoint_stem TEXT,
    config_path TEXT,
    index_path TEXT,
    embedding_backend TEXT,
    l3_embedding_path TEXT,
    metal_embedding_path TEXT,
    history_path TEXT,
    best_epoch INTEGER,
    val_mrr REAL,
    val_recall_at_1 REAL,
    val_recall_at_5 REAL,
    ga_set_id TEXT,
    ga_version_id TEXT,
    ga_binding_checksum TEXT,
    size_bytes INTEGER,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE,
    UNIQUE (project_id, run_id, checkpoint_path)
);

CREATE TABLE IF NOT EXISTS model_after_tasks (
    task_id TEXT PRIMARY KEY,
    task_dir TEXT NOT NULL,
    metal TEXT,
    embedding TEXT,
    candidate_file TEXT,
    positive_file TEXT,
    negative_file TEXT,
    notes TEXT,
    candidate_count INTEGER,
    positive_count INTEGER,
    negative_count INTEGER,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS model_after_batches (
    task_id TEXT NOT NULL REFERENCES model_after_tasks(task_id) ON DELETE CASCADE,
    batch_id TEXT NOT NULL,
    output_dir TEXT NOT NULL,
    command TEXT,
    model_runs_file TEXT,
    finished_at TEXT,
    n_models INTEGER,
    n_success INTEGER,
    best_model_id TEXT,
    summary_path TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (task_id, batch_id)
);

CREATE TABLE IF NOT EXISTS model_after_model_results (
    task_id TEXT NOT NULL,
    batch_id TEXT NOT NULL,
    model_id TEXT NOT NULL,
    model_run_root TEXT,
    checkpoint_path TEXT,
    project_id TEXT,
    run_id TEXT,
    registry_model_id TEXT,
    model_name TEXT,
    mrr REAL,
    hit_at_5 REAL,
    hit_at_10 REAL,
    hit_at_20 REAL,
    hit_at_50 REAL,
    selection_score REAL,
    rank_among_models INTEGER,
    status TEXT,
    ranking_path TEXT,
    metrics_path TEXT,
    report_path TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (task_id, batch_id, model_id),
    FOREIGN KEY (task_id, batch_id) REFERENCES model_after_batches(task_id, batch_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS jobs (
    job_id TEXT PRIMARY KEY,
    job_type TEXT NOT NULL,
    status TEXT NOT NULL,
    title TEXT,
    command_json TEXT,
    workspace_root TEXT,
    db_path TEXT,
    project_id TEXT,
    run_id TEXT,
    task_id TEXT,
    batch_id TEXT,
    model_id TEXT,
    log_path TEXT,
    result_json TEXT,
    error TEXT,
    created_at TEXT NOT NULL,
    started_at TEXT,
    finished_at TEXT,
    CHECK (status IN ('pending', 'running', 'success', 'failed', 'cancelled'))
);

CREATE TABLE IF NOT EXISTS hyperparameter_sets (
    hp_set_id TEXT PRIMARY KEY,
    hp_set_root TEXT NOT NULL,
    name TEXT,
    description TEXT,
    created_at TEXT,
    updated_at TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT
);

CREATE TABLE IF NOT EXISTS hyperparameter_versions (
    hp_set_id TEXT NOT NULL,
    hp_version_id TEXT NOT NULL,
    created_at TEXT,
    created_by TEXT,
    notes TEXT,
    yaml_path TEXT,
    source_json TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT,
    PRIMARY KEY (hp_set_id, hp_version_id),
    FOREIGN KEY (hp_set_id) REFERENCES hyperparameter_sets(hp_set_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_runs_project ON runs(project_id);
CREATE INDEX IF NOT EXISTS idx_step_executions_run ON step_executions(project_id, run_id);
CREATE INDEX IF NOT EXISTS idx_artifacts_run ON artifacts(project_id, run_id);
CREATE INDEX IF NOT EXISTS idx_models_run ON models(project_id, run_id);
CREATE INDEX IF NOT EXISTS idx_model_results_task ON model_after_model_results(task_id, batch_id);
CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON jobs(created_at);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_hp_versions_set ON hyperparameter_versions(hp_set_id, hp_version_id);
