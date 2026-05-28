"""Step definitions for Run-based ChemDB pipeline (Phase 1A core + Phase 1C training)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Sequence

ScriptKind = Literal["src_script", "repo_script", "python_module"]
RunnerKind = Literal["subprocess", "in_process"]


@dataclass(frozen=True)
class StepSpec:
    step_id: str
    script: str  # basename under src, repo-relative path, or module name (training.data)
    command_argv: List[str]  # template expanded by orchestrator; use placeholders
    required_inputs: List[str]  # relative to run_root (tmp/, data/, training/)
    expected_outputs: List[str]
    optional_outputs: List[str] = field(default_factory=list)
    log_name: Optional[str] = None  # defaults to step_id
    script_kind: ScriptKind = "src_script"
    runner: RunnerKind = "subprocess"

    @property
    def log_file(self) -> str:
        return f"logs/{self.log_name or self.step_id}.log"


PIPELINE_ORDER: List[str] = [
    "step1",
    "step2",
    "require_bound_ga",
    "apply_ga_to_run",
    "step6_7",
    "step8",
    "step12",
    "step13",
]

STEP_REGISTRY: Dict[str, StepSpec] = {
    "step1": StepSpec(
        step_id="step1",
        script="step1.py",
        command_argv=[
            "step1.py",
            "--pubchem-dir", "{data_dir}/pubchem",
            "--metal-list", "{data_dir}/metal_list.txt",
            "--tmp-dir", "{tmp_dir}",
        ],
        required_inputs=[
            "data/pubchem",
            "data/metal_list.txt",
        ],
        expected_outputs=[
            "tmp/complex_data.csv",
            "tmp/ligand_data.csv",
        ],
        optional_outputs=["tmp/processing_stats.json", "tmp/pubchem_data.csv"],
    ),
    "step2": StepSpec(
        step_id="step2",
        script="step2.py",
        command_argv=[
            "step2.py",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
        ],
        required_inputs=["tmp/ligand_data.csv"],
        expected_outputs=["tmp/repaired_ligand_data.csv"],
        optional_outputs=["tmp/repair_stats.json"],
    ),
    "require_bound_ga": StepSpec(
        step_id="require_bound_ga",
        script="",
        runner="in_process",
        command_argv=[],
        required_inputs=[
            "manifests/ga_binding.json",
            "tmp/GA_with_id.csv",
            "tmp/repaired_ligand_data.csv",
        ],
        expected_outputs=[],
    ),
    "apply_ga_to_run": StepSpec(
        step_id="apply_ga_to_run",
        script="step4_5.py",
        command_argv=[
            "step4_5.py",
            "--mode", "apply-ga",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
        ],
        required_inputs=[
            "tmp/repaired_ligand_data.csv",
            "tmp/GA_with_id.csv",
            "manifests/ga_binding.json",
        ],
        expected_outputs=[
            "tmp/ligand_with_gac.csv",
            "tmp/IRL_filtered.csv",
        ],
        optional_outputs=[
            "tmp/ligand_data_deduplicated.csv",
            "tmp/IRL_filtered_cleaned.csv",
            "tmp/step4_5_stats.json",
        ],
    ),
    "step4_5": StepSpec(
        step_id="step4_5",
        script="step4_5.py",
        command_argv=[
            "step4_5.py",
            "--mode", "full-auto",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
        ],
        required_inputs=["tmp/repaired_ligand_data.csv"],
        expected_outputs=[
            "tmp/GA_with_id.csv",
            "tmp/ligand_with_gac.csv",
            "tmp/IRL_filtered.csv",
        ],
        optional_outputs=[
            "tmp/ligand_data_deduplicated.csv",
            "tmp/IRL_filtered_cleaned.csv",
            "tmp/step4_5_stats.json",
        ],
        log_name="step4_5_legacy",
    ),
    "step6_7": StepSpec(
        step_id="step6_7",
        script="step6_7.py",
        command_argv=[
            "step6_7.py",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
        ],
        required_inputs=[
            "tmp/repaired_ligand_data.csv",
            "tmp/GA_with_id.csv",
            "tmp/IRL_filtered.csv",
        ],
        expected_outputs=["tmp/fragments.csv"],
        optional_outputs=["tmp/step6_7_stats.json"],
    ),
    "step8": StepSpec(
        step_id="step8",
        script="step8.py",
        command_argv=[
            "step8.py",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
            "--metal-list", "{data_dir}/metal_list.txt",
        ],
        required_inputs=[
            "tmp/complex_data.csv",
            "tmp/ligand_data.csv",
            "tmp/repaired_ligand_data.csv",
            "tmp/fragments.csv",
            "tmp/IRL_filtered.csv",
            "data/metal_list.txt",
        ],
        expected_outputs=[
            "tmp/neo4j_metals.csv",
            "tmp/neo4j_repaired_ligands.csv",
            "tmp/neo4j_l3_l4_relationships.csv",
            "tmp/neo4j_l4_l5_relationships.csv",
            "tmp/neo4j_m_l1_relationships.csv",
            "tmp/neo4j_l1_l3_relationships.csv",
        ],
        optional_outputs=[
            "tmp/neo4j_complexes.csv",
            "tmp/neo4j_ligands.csv",
            "tmp/neo4j_fragments.csv",
            "tmp/neo4j_irl.csv",
            "tmp/neo4j_l1_l2_relationships.csv",
            "tmp/neo4j_l2_l3_relationships.csv",
            "tmp/step8_stats.json",
        ],
    ),
    "step12": StepSpec(
        step_id="step12",
        script="step12.py",
        command_argv=[
            "step12.py",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
            "--K", "30",
        ],
        required_inputs=[
            "tmp/neo4j_l3_l4_relationships.csv",
            "tmp/neo4j_l4_l5_relationships.csv",
            "tmp/ligand_with_gac.csv",
            "tmp/neo4j_m_l1_relationships.csv",
            "tmp/neo4j_l1_l3_relationships.csv",
            "tmp/neo4j_repaired_ligands.csv",
        ],
        expected_outputs=[
            "tmp/l3_l5.json",
            "tmp/l5_l3.json",
            "tmp/l5_freq_weight.json",
            "tmp/l3_gac.json",
            "tmp/l5_l3_filtered_K30.json",
            "tmp/m_l3_pairs.csv",
        ],
        optional_outputs=["tmp/step12_stats.json"],
    ),
    "step13": StepSpec(
        step_id="step13",
        script="step_13_complete.py",
        command_argv=[
            "step_13_complete.py",
            "--input-dir", "{tmp_dir}",
            "--output-dir", "{tmp_dir}",
            "--kl-nl-only",
        ],
        required_inputs=[
            "tmp/m_l3_pairs.csv",
            "tmp/l3_l5.json",
            "tmp/l5_l3.json",
            "tmp/l5_freq_weight.json",
            "tmp/l3_gac.json",
        ],
        expected_outputs=[
            "tmp/step13_kl_nl_samples.csv",
        ],
        optional_outputs=[
            "tmp/step13_kl_nl_samples.stats.csv",
            "tmp/step13_stats.json",
            "tmp/metal_l3_index.csv",
        ],
    ),
}

TRAINING_PIPELINE_ORDER: List[str] = [
    "build_l3_embedding_ecfp",
    "zscore_metal_embedding",
    "prepare_training_index",
    "prepare_training_config",
    "training_data",
    "training_train",
]

TRAINING_STEP_REGISTRY: Dict[str, StepSpec] = {
    "build_l3_embedding_ecfp": StepSpec(
        step_id="build_l3_embedding_ecfp",
        script="tools/build_L3_embedding_index.py",
        command_argv=[
            "build_L3_embedding_index.py",
            "--tmp-dir", "{tmp_dir}",
            "--out-dir", "{data_dir}/L3_embedding",
            "--backends", "ecfp",
        ],
        required_inputs=[
            "tmp/repaired_ligand_data.csv",
            "tmp/metal_l3_index.csv",
        ],
        expected_outputs=["data/L3_embedding/L3_embedding_ecfp.npz"],
        optional_outputs=["data/L3_embedding/build.log"],
    ),
    "zscore_metal_embedding": StepSpec(
        step_id="zscore_metal_embedding",
        script="data/metal_embedding/zscore_element_features.py",
        script_kind="repo_script",
        command_argv=[
            "zscore_element_features.py",
            "--input", "{data_dir}/metal_embedding/element_features.csv",
            "--output", "{data_dir}/metal_embedding/element_features_zscore.csv",
        ],
        required_inputs=["data/metal_embedding/element_features.csv"],
        expected_outputs=["data/metal_embedding/element_features_zscore.csv"],
    ),
    "prepare_training_index": StepSpec(
        step_id="prepare_training_index",
        script="",
        runner="in_process",
        command_argv=[],
        required_inputs=[
            "tmp/step13_kl_nl_samples.csv",
            "tmp/m_l3_pairs.csv",
            "tmp/l3_gac.json",
            "data/L3_embedding/L3_embedding_ecfp.npz",
            "data/metal_embedding/element_features_zscore.csv",
        ],
        expected_outputs=["training/index.json"],
    ),
    "prepare_training_config": StepSpec(
        step_id="prepare_training_config",
        script="",
        runner="in_process",
        command_argv=[],
        required_inputs=["training/index.json"],
        expected_outputs=["training/config.yaml"],
    ),
    "training_data": StepSpec(
        step_id="training_data",
        script="training.data",
        script_kind="python_module",
        command_argv=[
            "training.data",
            "--training-dir", "{training_dir}",
            "--index", "{training_dir}/index.json",
            "--seed", "42",
        ],
        required_inputs=[
            "training/index.json",
            "tmp/step13_kl_nl_samples.csv",
            "tmp/m_l3_pairs.csv",
            "data/L3_embedding/L3_embedding_ecfp.npz",
            "data/metal_embedding/element_features_zscore.csv",
        ],
        expected_outputs=[
            "training/split_index.json",
            "training/train_records.pkl",
            "training/val_records.pkl",
            "training/test_records.pkl",
        ],
    ),
    "training_train": StepSpec(
        step_id="training_train",
        script="training.train",
        script_kind="python_module",
        command_argv=[
            "training.train",
            "--config", "{training_dir}/config.yaml",
        ],
        required_inputs=[
            "training/config.yaml",
            "training/train_records.pkl",
            "training/val_records.pkl",
            "training/test_records.pkl",
        ],
        expected_outputs=[
            "training/ckpts/best.pt",
            "training/ckpts/history.json",
        ],
        optional_outputs=["training/ckpts/last.pt"],
    ),
}

STEP_REGISTRY.update(TRAINING_STEP_REGISTRY)


def get_step(step_id: str) -> StepSpec:
    if step_id not in STEP_REGISTRY:
        raise KeyError(f"Unknown step_id: {step_id}")
    return STEP_REGISTRY[step_id]


def pipeline_through(through_step: str, *, pipeline: str = "core") -> List[StepSpec]:
    order = TRAINING_PIPELINE_ORDER if pipeline == "training" else PIPELINE_ORDER
    if through_step not in order:
        raise KeyError(f"through_step must be one of {order}, got {through_step}")
    idx = order.index(through_step) + 1
    return [STEP_REGISTRY[s] for s in order[:idx]]


def all_registered_paths_under_run(spec: StepSpec) -> List[str]:
    return list(spec.required_inputs) + list(spec.expected_outputs) + list(spec.optional_outputs)
