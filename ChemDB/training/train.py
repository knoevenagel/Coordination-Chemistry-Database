# -*- coding: utf-8 -*-
"""
训练入口：从统一配置文件加载超参，构建 data / model / loss，运行训练与验证，
按 val 指标保存 checkpoint，支持学习率调度与早停。
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import torch
import yaml

from . import data as data_module
from . import model as model_module

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 配置加载与默认值
# ---------------------------------------------------------------------------


def load_config(path: str | Path | None) -> dict:
    if path is None:
        base = Path(__file__).resolve().parent
        path = base / "config.yaml"
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not cfg:
        raise ValueError("Config is empty")
    return cfg


def resolve_config(cfg: dict) -> dict:
    base = Path(__file__).resolve().parent
    # data
    if cfg.get("data") is None:
        cfg["data"] = {}
    d = cfg["data"]
    if d.get("training_dir") is None:
        d["training_dir"] = str(base)
    if d.get("index_path") is None:
        d["index_path"] = str(base / "index.json")
    if d.get("embedding") is None:
        d["embedding"] = "gin"
    if d.get("batch_size") is None:
        d["batch_size"] = 32
    if d.get("num_workers") is None:
        d["num_workers"] = 0
    if d.get("seed") is None:
        d["seed"] = 42
    # model
    if cfg.get("model") is None:
        cfg["model"] = {}
    if cfg["model"].get("activation") is None:
        cfg["model"]["activation"] = "gelu"
    if cfg["model"].get("dropout") is None:
        cfg["model"]["dropout"] = 0.1
    # train
    if cfg.get("train") is None:
        cfg["train"] = {}
    t = cfg["train"]
    if t.get("epochs") is None:
        t["epochs"] = 50
    if t.get("lr") is None:
        t["lr"] = 1e-3
    if t.get("weight_decay") is None:
        t["weight_decay"] = 0.0
    if t.get("optimizer") is None:
        t["optimizer"] = "adamw"
    if t.get("scheduler") is None:
        t["scheduler"] = None
    if t.get("scheduler_patience") is None:
        t["scheduler_patience"] = 3
    if t.get("scheduler_factor") is None:
        t["scheduler_factor"] = 0.5
    if t.get("ckpt_dir") is None:
        t["ckpt_dir"] = str(base / "ckpts")
    if t.get("ckpt_best_metric") is None:
        t["ckpt_best_metric"] = "loss"
    if t.get("early_stop_patience") is None:
        t["early_stop_patience"] = 8
    if t.get("val_every") is None:
        t["val_every"] = 1
    if t.get("device") is None:
        t["device"] = "cuda"
    return cfg


def resolve_training_device(device_cfg: str | None) -> torch.device:
    """Use GPU for training; fail if CUDA is unavailable."""
    import os

    name = (device_cfg or os.environ.get("CHEMDB_TRAINING_DEVICE") or "cuda").strip()
    if name.lower() == "cpu":
        raise ValueError(
            "Training device must be GPU (cuda). Set train.device: cuda in config "
            "or CHEMDB_TRAINING_DEVICE=cuda:0."
        )
    if "cuda" in name.lower() and not torch.cuda.is_available():
        raise RuntimeError(
            f"Training requires GPU (device={name}) but torch.cuda.is_available() is False."
        )
    return torch.device(name)


# ---------------------------------------------------------------------------
# 验证指标：Recall@1, Recall@5, MRR
# ---------------------------------------------------------------------------


def compute_ranking_metrics(logits: torch.Tensor, target_idx: torch.Tensor) -> dict:
    """
    logits: (B, K), target_idx: (B,) 通常全 0。
    返回 recall_at_1, recall_at_5, mrr（标量）。
    """
    B = logits.shape[0]
    pred_rank = logits.argsort(dim=1, descending=True)
    rank_of_target = (pred_rank == target_idx.unsqueeze(1)).nonzero(as_tuple=True)
    if rank_of_target[0].numel() == 0:
        rank = torch.full((B,), logits.shape[1], device=logits.device, dtype=torch.long)
    else:
        rank = torch.full((B,), logits.shape[1], device=logits.device, dtype=torch.long)
        rank[rank_of_target[0]] = rank_of_target[1]
    recall_at_1 = (rank == 0).float().mean().item()
    recall_at_5 = (rank < 5).float().mean().item()
    mrr = (1.0 / (rank.float() + 1)).mean().item()
    return {"recall_at_1": recall_at_1, "recall_at_5": recall_at_5, "mrr": mrr}


# ---------------------------------------------------------------------------
# 训练与验证
# ---------------------------------------------------------------------------


def run_validation(model: torch.nn.Module, loader, device: torch.device, loss_fn):
    model.eval()
    total_loss = 0.0
    n = 0
    acc_r1 = 0.0
    acc_r5 = 0.0
    acc_mrr = 0.0
    with torch.no_grad():
        for batch in loader:
            metal = batch["metal"].to(device)
            context = batch["context"].to(device)
            candidates = batch["candidates"].to(device)
            target_idx = batch["target_idx"].to(device)
            logits = model(metal, context, candidates)
            loss = loss_fn(logits, target_idx)
            total_loss += loss.item() * metal.shape[0]
            n += metal.shape[0]
            m = compute_ranking_metrics(logits, target_idx)
            acc_r1 += m["recall_at_1"] * metal.shape[0]
            acc_r5 += m["recall_at_5"] * metal.shape[0]
            acc_mrr += m["mrr"] * metal.shape[0]
    if n == 0:
        return {"loss": 0.0, "recall_at_1": 0.0, "recall_at_5": 0.0, "mrr": 0.0}
    return {
        "loss": total_loss / n,
        "recall_at_1": acc_r1 / n,
        "recall_at_5": acc_r5 / n,
        "mrr": acc_mrr / n,
    }


def main():
    parser = argparse.ArgumentParser(description="条件排序模型训练（配置文件驱动）")
    parser.add_argument("--config", type=str, default=None, help="config.yaml 路径，默认 training/config.yaml")
    parser.add_argument("--resume", type=str, default=None, help="从 checkpoint 恢复（.pt 路径）")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    cfg = load_config(args.config)
    cfg = resolve_config(cfg)
    d_cfg = cfg["data"]
    m_cfg = cfg["model"]
    t_cfg = cfg["train"]

    device = resolve_training_device(t_cfg.get("device"))
    torch.manual_seed(d_cfg["seed"])

    paths = data_module.load_index_paths(d_cfg["index_path"])
    emb_key = f"ligand_embedding_{d_cfg['embedding']}"
    if emb_key not in paths:
        raise KeyError(f"Embedding {d_cfg['embedding']} not in index; keys: {list(paths.keys())}")
    ligand_path = paths[emb_key]
    metal_path = paths["metal_embedding"]

    train_loader = data_module.get_dataloader(
        split="train",
        training_dir=d_cfg["training_dir"],
        ligand_embedding_path=ligand_path,
        metal_embedding_path=metal_path,
        batch_size=d_cfg["batch_size"],
        max_candidates=d_cfg.get("max_candidates"),
        num_workers=d_cfg["num_workers"],
    )
    val_loader = data_module.get_dataloader(
        split="val",
        training_dir=d_cfg["training_dir"],
        ligand_embedding_path=ligand_path,
        metal_embedding_path=metal_path,
        batch_size=d_cfg["batch_size"],
        max_candidates=d_cfg.get("max_candidates"),
        num_workers=d_cfg["num_workers"],
    )

    ds = train_loader.dataset
    d_m = m_cfg.get("d_m")
    d_l = m_cfg.get("d_l")
    if d_m is None:
        d_m = ds.d_m
    if d_l is None:
        d_l = ds.d_l

    model = model_module.RankModel(
        d_m=d_m,
        d_l=d_l,
        activation=m_cfg["activation"],
        dropout=m_cfg["dropout"],
    ).to(device)

    loss_fn = torch.nn.CrossEntropyLoss()
    opt_name = (t_cfg.get("optimizer") or "adamw").lower()
    if opt_name == "adamw":
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=t_cfg["lr"],
            weight_decay=t_cfg["weight_decay"],
        )
    else:
        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=t_cfg["lr"],
            weight_decay=t_cfg["weight_decay"],
        )

    sched_cfg = t_cfg.get("scheduler")
    scheduler = None
    best_metric_key = t_cfg.get("ckpt_best_metric") or "loss"
    higher_is_better = best_metric_key not in ("loss", "val_loss")
    if sched_cfg == "reduce_on_plateau":
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode="max" if higher_is_better else "min",
            factor=t_cfg.get("scheduler_factor", 0.5),
            patience=t_cfg.get("scheduler_patience", 3),
        )
    elif sched_cfg == "cosine":
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=t_cfg["epochs"])

    ckpt_dir = Path(t_cfg["ckpt_dir"])
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    early_patience = int(t_cfg.get("early_stop_patience") or 0)
    val_every = int(t_cfg.get("val_every") or 1)

    start_epoch = 0
    best_metric_value = float("inf") if not higher_is_better else -1.0
    epochs_no_improve = 0
    history: list[dict] = []  # 每轮验证后追加：epoch, train_loss, val_loss, recall_at_1, recall_at_5, mrr

    if args.resume:
        ckpt = torch.load(args.resume, map_location=device)
        if "model_state_dict" in ckpt:
            model.load_state_dict(ckpt["model_state_dict"])
        else:
            model.load_state_dict(ckpt)
        start_epoch = ckpt.get("epoch", 0) + 1
        best_metric_value = ckpt.get("best_metric")
        if best_metric_value is None:
            best_metric_value = float("inf") if not higher_is_better else -1.0
        logger.info("Resumed from %s at epoch %s", args.resume, start_epoch)

    for epoch in range(start_epoch, t_cfg["epochs"]):
        model.train()
        train_loss = 0.0
        train_n = 0
        for batch in train_loader:
            metal = batch["metal"].to(device)
            context = batch["context"].to(device)
            candidates = batch["candidates"].to(device)
            target_idx = batch["target_idx"].to(device)
            optimizer.zero_grad()
            logits = model(metal, context, candidates)
            loss = loss_fn(logits, target_idx)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * metal.shape[0]
            train_n += metal.shape[0]
        train_loss = train_loss / train_n if train_n else 0.0

        if (epoch + 1) % val_every == 0 or epoch == start_epoch:
            val_metrics = run_validation(model, val_loader, device, loss_fn)
            current_best = val_metrics.get(best_metric_key, val_metrics.get("loss", val_metrics["mrr"]))
            if scheduler and sched_cfg == "reduce_on_plateau":
                scheduler.step(current_best)
            elif scheduler and sched_cfg == "cosine":
                scheduler.step()

            logger.info(
                "epoch %d train_loss=%.4f val_loss=%.4f val_R@1=%.4f val_R@5=%.4f val_MRR=%.4f",
                epoch + 1,
                train_loss,
                val_metrics["loss"],
                val_metrics["recall_at_1"],
                val_metrics["recall_at_5"],
                val_metrics["mrr"],
            )
            history.append({
                "epoch": epoch + 1,
                "train_loss": float(train_loss),
                "val_loss": float(val_metrics["loss"]),
                "recall_at_1": float(val_metrics["recall_at_1"]),
                "recall_at_5": float(val_metrics["recall_at_5"]),
                "mrr": float(val_metrics["mrr"]),
            })

            improved = current_best < best_metric_value if not higher_is_better else current_best > best_metric_value
            if improved:
                best_metric_value = current_best
                epochs_no_improve = 0
                ckpt_path = ckpt_dir / "best.pt"
                torch.save(
                    {
                        "model_state_dict": model.state_dict(),
                        "epoch": epoch,
                        "best_metric": best_metric_value,
                        "config": cfg,
                    },
                    ckpt_path,
                )
                logger.info("Saved best checkpoint to %s (%s=%.4f)", ckpt_path, best_metric_key, best_metric_value)
            else:
                epochs_no_improve += 1

            if early_patience > 0 and epochs_no_improve >= early_patience:
                logger.info("Early stopping at epoch %d (no improvement for %d validations)", epoch + 1, early_patience)
                break
        else:
            if scheduler and sched_cfg == "cosine":
                scheduler.step()

    # last.pt 使用实际跑到的最后一轮（含早停情况）
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "epoch": epoch,
            "config": cfg,
        },
        ckpt_dir / "last.pt",
    )
    # 保存 epoch 与性能关系，便于画曲线或对比
    history_path = ckpt_dir / "history.json"
    with open(history_path, "w", encoding="utf-8") as f:
        json.dump(history, f, indent=2, ensure_ascii=False)
    logger.info("Training done. Last checkpoint: %s, history: %s", ckpt_dir / "last.pt", history_path)


if __name__ == "__main__":
    main()
    sys.exit(0)
