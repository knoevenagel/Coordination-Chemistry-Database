"""Run-based model-after: load bundle, rank candidates, evaluate task fit, multi-model selection."""

from .bundle import RunModelBundle, load_run_bundle
from .batch import evaluate_models
from .eval_task import evaluate_task_fit
from .io import TaskBundle, load_task_bundle
from .scoring import recommend_single_model

__all__ = [
    "RunModelBundle",
    "load_run_bundle",
    "TaskBundle",
    "load_task_bundle",
    "recommend_single_model",
    "evaluate_task_fit",
    "evaluate_models",
]
