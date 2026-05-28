#!/usr/bin/env python3
"""
ChemDB Interactive Entry Point
"""

import os
import sys
import time
import logging
import traceback
import humanize
from pathlib import Path
from datetime import datetime
from multiprocessing import cpu_count

from InquirerPy import inquirer
from InquirerPy.separator import Separator
from InquirerPy.validator import NumberValidator
from InquirerPy.utils import color_print

# Add src to path for imports
_script_dir = Path(__file__).parent if "__file__" in dir() else Path.cwd()
sys.path.insert(0, str(_script_dir / "src"))

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Step definitions with output files for "last generated" detection
STEP_DEFINITIONS = {
    "step1": {
        "name": "Step 1: Initial data processing",
        "output_files": ["./tmp/complex_data.csv", "./tmp/ligand_data.csv"],
        "processor_class": "Step1Processor",
        "module": "step1",
        "method": "process_all_files",
        "supports_workers": True,
        "supports_limit": True,
    },
    "step2": {
        "name": "Step 2: Ligand repair",
        "output_files": ["./tmp/repaired_ligand_data.csv"],
        "processor_class": "Step2Processor",
        "module": "step2",
        "method": "process_all_ligands",
        "supports_workers": True,
        "supports_limit": True,
    },
    "step4_5": {
        "name": "Step 4/5: Generate IRL and GA",
        "output_files": ["./tmp/IRL_filtered.csv", "./tmp/GA_with_id.csv"],
        "processor_class": "Step4_5Processor",
        "module": "step4_5",
        "method": "process_all",
        "supports_workers": True,
        "supports_limit": True,
    },
    "step6_7": {
        "name": "Step 6/7: Ligand marking and L4 fragment generation",
        "output_files": ["./tmp/fragments.csv"],
        "processor_class": "Step6_7Processor",
        "module": "step6_7",
        "method": "process_all_ligands",
        "supports_workers": True,
        "supports_limit": True,
    },
    "step8": {
        "name": "Step 8: Generate CSV files for Neo4j bulk import",
        "output_files": ["./tmp/neo4j_metals.csv", "./tmp/neo4j_complexes.csv"],
        "processor_class": "Step8Processor",
        "module": "step8",
        "method": "process_all_data",
        "supports_workers": False,
        "supports_limit": True,
    },
    "step9": {
        "name": "Step 9: Generate SMILES fingerprint",
        "output_files": ["./tmp/did_index.csv"],
        "processor_class": "Step9Processor",
        "module": "step9",
        "method": "process_all_smiles",
        "supports_workers": True,
        "supports_limit": True,
    },
    # "step10": {
    #     "name": "Step 10: Merge layers and generate combined database",
    #     "output_files": ["./tmp/step10_layers.csv"],
    #     "processor_class": "Step10Processor",
    #     "module": "step10",
    #     "method": "process",
    #     "supports_workers": False,
    #     "supports_limit": True,
    # },
    "step11": {
        "name": "Step 11: Generate embeddings of L3 and M-L3 relationships",
        "output_files": ["./tmp/l3_embeddings.csv", "./tmp/m_l3_relationships.csv"],
        "processor_class": "Step11Processor",
        "module": "step11",
        "method": "process_all",
        "supports_workers": False,
        "supports_limit": True,
        "supports_model_type": True,
    },
}

# Ordered list of step keys
STEP_ORDER = ["step1", "step2", "step4_5", "step6_7", "step8", "step9", "step11"]


def get_last_modified_time(file_path: str) -> str:
    """Get the last modification time of a file as a formatted string."""
    path = Path(file_path)
    if path.exists():
        mtime = path.stat().st_mtime
        dt = datetime.fromtimestamp(mtime)
        return humanize.naturaltime(dt)
        # return dt.strftime("%Y-%m-%d %H:%M:%S")
    return None


def get_step_output_status(step_key: str) -> str:
    """Check if step output exists and return status string."""
    step_def = STEP_DEFINITIONS[step_key]
    output_files = step_def["output_files"]
    
    # Check which files exist and find the most recent modification time
    existing_files = []
    latest_time = None
    
    for output_file in output_files:
        time_str = get_last_modified_time(output_file)
        if time_str:
            existing_files.append(output_file)
            if latest_time is None:
                latest_time = time_str
    
    # Return appropriate status
    if not existing_files:
        return "none"  # All files missing
    elif len(existing_files) < len(output_files):
        return "missing"  # Some files missing
    else:
        return latest_time  # All files exist


def build_step_choices():
    """Build choices for step selection with output status."""
    choices = []
    for step_key in STEP_ORDER:
        step_def = STEP_DEFINITIONS[step_key]
        # status = get_step_output_status(step_key)
        name = step_def['name']
        choices.append({"name": name, "value": step_key})
    return choices


def run_step(step_key: str, num_workers: int, record_limit: int, model_type: str = "gin") -> bool:
    """Run a single step and return success status."""
    step_def = STEP_DEFINITIONS[step_key]
    
    logger.info("=" * 60)
    logger.info(f"STARTING {step_def['name'].upper()}")
    logger.info("=" * 60)
    
    try:
        # Import the module dynamically
        module = __import__(step_def["module"])
        processor_class = getattr(module, step_def["processor_class"])
        
        # Build constructor arguments based on step support
        kwargs = {}
        if step_def.get("supports_workers", False) and num_workers is not None:
            kwargs["num_workers"] = num_workers
        if step_def.get("supports_limit", False):
            kwargs["record_limit"] = record_limit
        if step_def.get("supports_model_type", False):
            kwargs["model_type"] = model_type
        
        # Create processor and run
        start_time = time.time()
        processor = processor_class(**kwargs)
        
        method = getattr(processor, step_def["method"])
        stats = method()
        
        elapsed = time.time() - start_time
        logger.info(f"{step_def['name']} completed successfully in {elapsed:.2f} seconds")
        
        return True
        
    except Exception as e:
        logger.error(f"{step_def['name']} failed: {e}")
        logger.error(f"Error details: {traceback.format_exc()}")
        return False


def run_selected_steps(selected_steps: list, num_workers: int, record_limit: int, model_type: str = "gin") -> bool:
    """Run selected steps in order, stopping on failure."""
    # Sort steps by their order in STEP_ORDER
    ordered_steps = [s for s in STEP_ORDER if s in selected_steps]
    
    logger.info("=" * 60)
    logger.info("STARTING SEQUENTIAL EXECUTION OF SELECTED STEPS")
    logger.info(f"Steps to run: {', '.join(ordered_steps)}")
    logger.info(f"Workers: {num_workers if num_workers else 'Auto'}")
    logger.info(f"Record Limit: {record_limit if record_limit > 0 else 'Unlimited'}")
    logger.info(f"Model Type (for step11): {model_type}")
    logger.info("=" * 60)
    
    overall_start = time.time()
    
    for step_key in ordered_steps:
        success = run_step(step_key, num_workers, record_limit, model_type)
        if not success:
            logger.error(f"{step_key} failed. Stopping execution.")
            return False
    
    overall_time = time.time() - overall_start
    
    logger.info("=" * 60)
    logger.info("ALL SELECTED STEPS COMPLETED SUCCESSFULLY!")
    logger.info(f"Total execution time: {overall_time:.2f} seconds")
    logger.info("=" * 60)
    
    return True


def main():
    """Main interactive entry point."""
    print("=" * 60)
    print("ChemDB Interactive Pipeline")
    print("=" * 60)
    
    # 1. Multi-select steps

    print("Steps status:")
    for step_key in STEP_DEFINITIONS.keys():
        status = get_step_output_status(step_key)
        if status == "none":
            color_print(formatted_text=[("class:aa", f"  {step_key}: "), ("class:bb", "No files generated")], style={"aa": "white", "bb": "grey"})
        elif status == "missing":
            color_print(formatted_text=[("class:aa", f"  {step_key}: "), ("class:bb", "Some files missing")], style={"aa": "white", "bb": "red"})
        else:
            color_print(formatted_text=[("class:aa", f"  {step_key}: "), ("class:bb", f"Last generated {status}")], style={"aa": "white", "bb": "green"})

    step_choices = build_step_choices()
    choices = [{"name": "[Select All]", "value": "__all__"}, Separator()] + step_choices

    selected_steps = inquirer.checkbox(
        message="Select steps to run:",
        choices=choices,
        instruction="space to select, enter to confirm",
        transformer=lambda result: f"{'All' if '[Select All]' in result else len(result)} step(s) selected",
    ).execute()
    
    if not selected_steps:
        print("No steps selected. Exiting.")
        return 0
    
    if "__all__" in selected_steps:
        selected_steps = [c["value"] for c in step_choices]
    
    # 2. Record limit input
    record_limit_str = inquirer.text(
        message="Enter record limit (leave empty for unlimited):",
        default="",
    ).execute()
    
    if record_limit_str.strip() == "":
        record_limit = -1
    else:
        try:
            record_limit = int(record_limit_str)
            if record_limit <= 0:
                record_limit = -1
        except ValueError:
            print("Invalid number, using unlimited.")
            record_limit = -1
    
    # 3. Worker number input
    default_workers = cpu_count()
    worker_str = inquirer.text(
        message=f"Enter number of workers (leave empty for auto={default_workers}):",
        default="",
    ).execute()
    
    if worker_str.strip() == "":
        num_workers = None  # Auto
    else:
        try:
            num_workers = int(worker_str)
            if num_workers <= 0:
                num_workers = None
        except ValueError:
            print("Invalid number, using auto.")
            num_workers = None
    
    # 4. Model type for step11 (only if step11 is selected)
    model_type = "gin"
    if "step11" in selected_steps:
        model_type = inquirer.select(
            message="Select model type for Step 11:",
            choices=[
                {"name": "GIN (Graph Isomorphism Network)", "value": "gin"},
                {"name": "GCN (Graph Convolutional Network)", "value": "gcn"},
            ],
            default="gin",
        ).execute()
    
    # 5. Confirmation
    print("\n" + "=" * 60)
    print("Configuration Summary:")
    print("=" * 60)
    print(f"Selected steps: {', '.join(selected_steps)}")
    print(f"Record limit: {record_limit if record_limit > 0 else 'Unlimited'}")
    print(f"Workers: {num_workers if num_workers else 'Auto'}")
    if "step11" in selected_steps:
        print(f"Model type: {model_type}")
    print("=" * 60)
    
    confirm = inquirer.confirm(
        message="Proceed with execution?",
        default=True
    ).execute()
    
    if not confirm:
        print("Execution cancelled.")
        return 0
    
    # 6. Run the pipeline
    try:
        success = run_selected_steps(selected_steps, num_workers, record_limit, model_type)
        return 0 if success else 1
    except KeyboardInterrupt:
        logger.info("Execution interrupted by user")
        print("\nExecution was interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Execution failed: {e}")
        logger.error(f"Error details: {traceback.format_exc()}")
        print(f"\nExecution failed: {e}")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
