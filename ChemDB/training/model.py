"""
条件排序打分网络：Query 编码器（MLP） + 点积打分。
结构固定：LayerNorm → 4×d_l → 2×d_l → d_l；logits = q · candidates。
"""
from __future__ import annotations

import torch
import torch.nn as nn
from typing import Optional

# 激活函数名 -> 模块
_ACT = {
    "gelu": nn.GELU,
    "relu": nn.ReLU,
    "silu": nn.SiLU,
}


def _get_activation(name: str):
    name = (name or "gelu").strip().lower()
    if name not in _ACT:
        raise ValueError(f"Unknown activation: {name}. Choose from {list(_ACT.keys())}")
    return _ACT[name]()


class RankModel(nn.Module):
    """
    输入: metal (B, d_m), context (B, d_l), candidates (B, K, d_l)
    输出: logits (B, K)，与 algorithm 中 f(M, K_L, l) 对应；CrossEntropy 目标为 target_idx=0。
    """

    def __init__(
        self,
        d_m: int,
        d_l: int,
        activation: str = "gelu",
        dropout: float = 0.1,
    ):
        super().__init__()
        self.d_m = d_m
        self.d_l = d_l
        d_in = d_m + d_l
        act_fn = _get_activation(activation)

        self.norm = nn.LayerNorm(d_in)
        self.fc1 = nn.Linear(d_in, 4 * d_l)
        self.act1 = act_fn
        self.drop1 = nn.Dropout(p=dropout)
        self.fc2 = nn.Linear(4 * d_l, 2 * d_l)
        self.act2 = act_fn
        self.drop2 = nn.Dropout(p=dropout)
        self.fc3 = nn.Linear(2 * d_l, d_l)

    def forward(
        self,
        metal: torch.Tensor,
        context: torch.Tensor,
        candidates: torch.Tensor,
    ) -> torch.Tensor:
        """
        metal: (B, d_m), context: (B, d_l), candidates: (B, K, d_l)
        returns: logits (B, K)
        """
        x = torch.cat([metal, context], dim=-1)
        x = self.norm(x)
        x = self.drop1(self.act1(self.fc1(x)))
        x = self.drop2(self.act2(self.fc2(x)))
        q = self.fc3(x)
        logits = torch.einsum("bd,bkd->bk", q, candidates)
        return logits
