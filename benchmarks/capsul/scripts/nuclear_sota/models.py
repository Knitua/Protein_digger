from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class ModelOutput:
    logits: torch.Tensor
    attention: torch.Tensor
    query_diversity_loss: torch.Tensor
    spectral_loss: torch.Tensor
    evidence: torch.Tensor | None = None
    auxiliary: dict[str, torch.Tensor] | None = None


def smooth_attention_logits(logits: torch.Tensor, radius: int) -> torch.Tensor:
    """Normalized local smoothing over residues, preserving masked logits downstream."""
    if radius <= 0:
        return logits
    shape = logits.shape
    flattened = logits.reshape(-1, 1, shape[-1])
    smoothed = F.avg_pool1d(
        flattened, kernel_size=2 * radius + 1, stride=1, padding=radius, count_include_pad=False
    )
    return smoothed.reshape(shape)


def dct_high_frequency_energy(attention: torch.Tensor, cutoff_fraction: float = 0.25) -> torch.Tensor:
    """Differentiable DCT-II high-frequency energy using an even FFT extension."""
    length = attention.shape[-1]
    if length < 8:
        return attention.new_zeros(())
    extension = torch.cat((attention, torch.flip(attention, dims=(-1,))), dim=-1)
    spectrum = torch.fft.rfft(extension.float(), dim=-1)[..., :length]
    frequency = torch.arange(length, device=attention.device, dtype=torch.float32)
    phase = torch.exp(-0.5j * math.pi * frequency / length)
    dct = (spectrum * phase).real
    start = max(1, int(math.ceil(length * cutoff_fraction)))
    numerator = dct[..., start:].square().mean()
    denominator = dct[..., 1:].square().mean().clamp_min(1e-12)
    return numerator / denominator


def query_diversity_loss(query: torch.Tensor) -> torch.Tensor:
    """Penalize off-diagonal cosine similarity between label queries."""
    if query.shape[0] <= 1:
        return query.new_zeros(())
    flattened = F.normalize(query.flatten(1), dim=-1)
    gram = flattened @ flattened.transpose(0, 1)
    off_diagonal = gram - torch.eye(gram.shape[0], device=gram.device, dtype=gram.dtype)
    return off_diagonal.square().sum() / (gram.numel() - gram.shape[0])


class LabelQueryHead(nn.Module):
    """Multi-head residue pooling with shared, label-specific, or dual label queries."""

    def __init__(
        self,
        input_dim: int,
        projection_dim: int,
        output_dim: int,
        heads: int = 4,
        mode: str = "label",
        dropout: float = 0.15,
        smoothing_radius: int = 0,
        dct_cutoff_fraction: float = 0.25,
    ) -> None:
        super().__init__()
        if projection_dim % heads:
            raise ValueError("projection_dim must be divisible by heads")
        if mode not in {"shared", "label", "label_local_global"}:
            raise ValueError(f"Unknown query mode: {mode}")
        self.output_dim = output_dim
        self.heads = heads
        self.head_dim = projection_dim // heads
        self.mode = mode
        self.smoothing_radius = smoothing_radius
        self.dct_cutoff_fraction = dct_cutoff_fraction
        query_labels = 1 if mode == "shared" else output_dim
        query_slots = 2 if mode == "label_local_global" else 1
        self.query = nn.Parameter(torch.empty(query_labels, query_slots, heads, self.head_dim))
        nn.init.normal_(self.query, std=1.0 / math.sqrt(self.head_dim))
        self.key = nn.Linear(input_dim, projection_dim, bias=False)
        self.value = nn.Linear(input_dim, projection_dim, bias=False)
        self.key_norm = nn.LayerNorm(self.head_dim)
        self.value_norm = nn.LayerNorm(self.head_dim)
        feature_dim = projection_dim * query_slots
        hidden_dim = max(128, projection_dim // 2)
        self.label_norm = nn.LayerNorm(feature_dim)
        self.label_mlp = nn.Sequential(
            nn.Linear(feature_dim, projection_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(projection_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.classifier_weight = nn.Parameter(torch.empty(output_dim, hidden_dim))
        self.classifier_bias = nn.Parameter(torch.zeros(output_dim))
        nn.init.xavier_uniform_(self.classifier_weight)

    def forward(
        self,
        residues: torch.Tensor,
        mask: torch.Tensor,
        attention_prior: torch.Tensor | None = None,
    ) -> ModelOutput:
        batch_size, length, _ = residues.shape
        key = self.key(residues).reshape(batch_size, length, self.heads, self.head_dim)
        value = self.value(residues).reshape(batch_size, length, self.heads, self.head_dim)
        key = self.key_norm(key)
        value = self.value_norm(value)
        query = self.query
        if query.shape[0] == 1:
            query = query.expand(self.output_dim, -1, -1, -1)
        scores = torch.einsum("blhd,cshd->bcshl", key, query) / math.sqrt(self.head_dim)
        scores = smooth_attention_logits(scores, self.smoothing_radius)
        if attention_prior is not None:
            if attention_prior.shape != (batch_size, self.output_dim, length):
                raise ValueError("attention_prior must have shape [batch, label, length]")
            scores = scores + attention_prior[:, :, None, None, :].to(scores.dtype)
        scores = scores.float().masked_fill(
            ~mask[:, None, None, None, :], torch.finfo(torch.float32).min
        )
        attention = torch.softmax(scores, dim=-1).to(value.dtype)
        if self.mode == "label_local_global":
            uniform = mask[:, None, None, :].to(attention.dtype)
            uniform = uniform / uniform.sum(dim=-1, keepdim=True).clamp_min(1.0)
            uniform = uniform.expand(-1, self.output_dim, self.heads, -1)
            attention = attention.clone()
            attention[:, :, 1] = 0.5 * attention[:, :, 1] + 0.5 * uniform
        pooled = torch.einsum("bcshl,blhd->bcshd", attention, value)
        pooled = pooled.reshape(batch_size, self.output_dim, -1).float()
        hidden = self.label_mlp(self.label_norm(pooled))
        logits = torch.einsum("bch,ch->bc", hidden, self.classifier_weight) + self.classifier_bias
        label_attention = attention.mean(dim=(2, 3))
        return ModelOutput(
            logits=logits,
            attention=label_attention,
            query_diversity_loss=query_diversity_loss(query),
            spectral_loss=dct_high_frequency_energy(label_attention, self.dct_cutoff_fraction),
            evidence=hidden,
            auxiliary=None,
        )


class LegacyN08Head(nn.Module):
    """Exact shared residue pooling architecture used by the frozen N08 baseline."""

    def __init__(
        self,
        input_dim: int,
        projection_dim: int,
        output_dim: int,
        dropout: float = 0.25,
        residual_blocks: int = 1,
    ) -> None:
        super().__init__()
        self.query = nn.Parameter(torch.empty(input_dim))
        nn.init.normal_(self.query, mean=0.0, std=1.0 / math.sqrt(input_dim))
        self.pool_norm = nn.LayerNorm(input_dim)
        self.in_proj = nn.Sequential(
            nn.Linear(input_dim, projection_dim), nn.GELU(), nn.Dropout(dropout)
        )
        self.blocks = nn.ModuleList(
            [
                nn.Sequential(
                    nn.LayerNorm(projection_dim),
                    nn.Linear(projection_dim, projection_dim * 2),
                    nn.GELU(),
                    nn.Dropout(dropout),
                    nn.Linear(projection_dim * 2, projection_dim),
                    nn.Dropout(dropout),
                )
                for _ in range(residual_blocks)
            ]
        )
        self.out = nn.Sequential(
            nn.LayerNorm(projection_dim),
            nn.Linear(projection_dim, max(64, projection_dim // 2)),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(max(64, projection_dim // 2), output_dim),
        )

    def forward(self, residues: torch.Tensor, mask: torch.Tensor) -> ModelOutput:
        scores = torch.einsum("bld,d->bl", residues, self.query) / math.sqrt(residues.shape[-1])
        scores = scores.float().masked_fill(~mask, torch.finfo(torch.float32).min)
        attention = torch.softmax(scores, dim=-1).to(residues.dtype)
        pooled = torch.einsum("bl,bld->bd", attention, residues).float()
        hidden = self.in_proj(self.pool_norm(pooled))
        for block in self.blocks:
            hidden = hidden + block(hidden)
        logits = self.out(hidden)
        label_attention = attention[:, None, :].expand(-1, logits.shape[1], -1)
        zero = logits.new_zeros(())
        return ModelOutput(logits, label_attention, zero, zero, None, None)


class NucleusSpecialistHead(nn.Module):
    """Three-query import/export/negative-evidence nucleus specialist."""

    def __init__(self, input_dim: int, projection_dim: int = 384, heads: int = 4, dropout: float = 0.15):
        super().__init__()
        if projection_dim % heads:
            raise ValueError("projection_dim must be divisible by heads")
        self.heads = heads
        self.head_dim = projection_dim // heads
        self.query = nn.Parameter(torch.empty(3, heads, self.head_dim))
        nn.init.normal_(self.query, std=1.0 / math.sqrt(self.head_dim))
        self.key = nn.Linear(input_dim, projection_dim, bias=False)
        self.value = nn.Linear(input_dim, projection_dim, bias=False)
        self.expert_score = nn.Sequential(
            nn.LayerNorm(projection_dim),
            nn.Linear(projection_dim, projection_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(projection_dim // 2, 1),
        )
        self.bias = nn.Parameter(torch.zeros(()))

    def forward(self, residues: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch, length, _ = residues.shape
        key = self.key(residues).reshape(batch, length, self.heads, self.head_dim)
        value = self.value(residues).reshape(batch, length, self.heads, self.head_dim)
        score = torch.einsum("blhd,qhd->bqhl", key, self.query) / math.sqrt(self.head_dim)
        score = score.float().masked_fill(~mask[:, None, None, :], torch.finfo(torch.float32).min)
        attention = torch.softmax(score, dim=-1).to(value.dtype)
        pooled = torch.einsum("bqhl,blhd->bqhd", attention, value).reshape(batch, 3, -1)
        expert = self.expert_score(pooled.float()).squeeze(-1)
        logit = expert[:, 0] + expert[:, 1] - expert[:, 2] + self.bias
        return logit, attention.mean(dim=2)


class GeneralWithNucleusSpecialist(nn.Module):
    def __init__(
        self,
        general_head: LabelQueryHead,
        input_dim: int,
        nucleus_index: int,
        specialist_projection_dim: int = 384,
        specialist_heads: int = 4,
        dropout: float = 0.15,
    ) -> None:
        super().__init__()
        self.general_head = general_head
        self.nucleus_index = nucleus_index
        self.specialist = NucleusSpecialistHead(
            input_dim, specialist_projection_dim, specialist_heads, dropout
        )
        self.fusion_logit = nn.Parameter(torch.tensor(-1.0))

    def forward(self, residues: torch.Tensor, mask: torch.Tensor) -> ModelOutput:
        general = self.general_head(residues, mask)
        specialist_logit, specialist_attention = self.specialist(residues, mask)
        generic_nucleus_logit = general.logits[:, self.nucleus_index]
        fusion_weight = torch.sigmoid(self.fusion_logit)
        logits = general.logits.clone()
        logits[:, self.nucleus_index] = generic_nucleus_logit + fusion_weight * specialist_logit
        spectral = general.spectral_loss + dct_high_frequency_energy(specialist_attention)
        return ModelOutput(
            logits=logits,
            attention=general.attention,
            query_diversity_loss=general.query_diversity_loss + query_diversity_loss(self.specialist.query),
            spectral_loss=spectral,
            evidence=general.evidence,
            auxiliary={
                "generic_nucleus_logit": generic_nucleus_logit,
                "specialist_nucleus_logit": specialist_logit,
                "specialist_attention": specialist_attention,
                "fusion_weight": fusion_weight,
            },
        )


class WindowMILHead(nn.Module):
    """Aggregate per-window label evidence with position-aware protein-level MIL."""

    def __init__(self, residue_head: LabelQueryHead, aggregation: str = "attention") -> None:
        super().__init__()
        if aggregation not in {"mean", "attention", "top2_noisy_or"}:
            raise ValueError(f"Unsupported MIL aggregation: {aggregation}")
        self.residue_head = residue_head
        self.aggregation = aggregation
        hidden_dim = residue_head.classifier_weight.shape[1]
        output_dim = residue_head.output_dim
        self.metadata_projection = nn.Sequential(
            nn.Linear(5, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim)
        )
        self.gate_weight = nn.Parameter(torch.empty(output_dim, hidden_dim))
        self.gate_bias = nn.Parameter(torch.zeros(output_dim))
        nn.init.xavier_uniform_(self.gate_weight)

    def classify(self, evidence: torch.Tensor) -> torch.Tensor:
        return (
            torch.einsum("bch,ch->bc", evidence, self.residue_head.classifier_weight)
            + self.residue_head.classifier_bias
        )

    def forward(
        self,
        windows: torch.Tensor,
        window_mask: torch.Tensor,
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
        window_valid: torch.Tensor | None = None,
    ) -> ModelOutput:
        if window_valid is None:
            window_valid = torch.ones(
                windows.shape[0], dtype=torch.bool, device=windows.device
            )
        if window_valid.shape != (windows.shape[0],):
            raise ValueError("window_valid must have shape [window]")
        window_output = self.residue_head(windows, window_mask)
        if window_output.evidence is None:
            raise RuntimeError("MIL requires per-label evidence tokens")
        evidence = window_output.evidence + self.metadata_projection(metadata.float())[:, None, :]
        protein_evidence = []
        protein_logits = []
        protein_availability = []
        window_gate = evidence.new_zeros((windows.shape[0], self.residue_head.output_dim))
        for protein_index in range(protein_count):
            protein_windows = owner == protein_index
            if not protein_windows.any():
                raise ValueError(f"Protein {protein_index} has no windows")
            selected = protein_windows & window_valid
            selected_evidence = evidence[selected]
            if selected_evidence.shape[0] == 0:
                # An expert may be entirely unavailable for one protein (for
                # example, no confident 3Di residues).  Emit a neutral token;
                # the downstream expert-availability mask excludes it.
                neutral_evidence = evidence[protein_windows][0] * 0.0
                protein_availability.append(False)
                if self.aggregation == "top2_noisy_or":
                    protein_logits.append(self.classify(neutral_evidence[None, ...])[0] * 0.0)
                else:
                    protein_evidence.append(neutral_evidence)
                continue
            protein_availability.append(True)
            if self.aggregation == "mean":
                aggregate = selected_evidence.mean(dim=0)
                window_gate[selected] = 1.0 / selected_evidence.shape[0]
                protein_evidence.append(aggregate)
            elif self.aggregation == "attention":
                gate_logits = torch.einsum(
                    "wch,ch->wc", selected_evidence, self.gate_weight
                ) + self.gate_bias
                gate = torch.softmax(gate_logits, dim=0)
                # Autocast can leave the softmax in fp32 while the evidence-backed
                # output buffer is fp16/bf16.  Keep the numerically safer softmax,
                # but cast at the two mixed-precision boundaries explicitly.
                gate_for_evidence = gate.to(selected_evidence.dtype)
                window_gate[selected] = gate.to(window_gate.dtype)
                aggregate = torch.einsum(
                    "wc,wch->ch", gate_for_evidence, selected_evidence
                )
                protein_evidence.append(aggregate)
            else:
                selected_probability = torch.sigmoid(self.classify(selected_evidence))
                top_count = min(2, selected_probability.shape[0])
                top = torch.topk(selected_probability, top_count, dim=0)
                top_probability = top.values
                selected_gate = torch.zeros_like(selected_probability)
                selected_gate.scatter_(0, top.indices, 1.0 / top_count)
                window_gate[selected] = selected_gate
                noisy_or = 1.0 - torch.prod(1.0 - top_probability, dim=0)
                noisy_or = noisy_or.clamp(1e-6, 1.0 - 1e-6)
                protein_logits.append(torch.logit(noisy_or))
        if self.aggregation == "top2_noisy_or":
            logits = torch.stack(protein_logits)
            aggregated_evidence = None
        else:
            aggregated_evidence = torch.stack(protein_evidence)
            logits = self.classify(aggregated_evidence)
        return ModelOutput(
            logits=logits,
            attention=window_output.attention,
            query_diversity_loss=window_output.query_diversity_loss,
            spectral_loss=window_output.spectral_loss,
            evidence=aggregated_evidence,
            auxiliary={
                "protein_availability": torch.tensor(
                    protein_availability, dtype=torch.bool, device=windows.device
                ),
                "window_gate": window_gate,
            },
        )


class WindowMILWithNucleusSpecialist(nn.Module):
    """Full-window general MIL plus an independently aggregated nucleus specialist."""

    def __init__(
        self,
        residue_head: LabelQueryHead,
        aggregation: str,
        nucleus_index: int,
        specialist_projection_dim: int = 384,
        specialist_heads: int = 4,
        dropout: float = 0.15,
    ) -> None:
        super().__init__()
        self.general_mil = WindowMILHead(residue_head, aggregation)
        self.specialist = NucleusSpecialistHead(
            residue_head.key.in_features,
            specialist_projection_dim,
            specialist_heads,
            dropout,
        )
        self.aggregation = aggregation
        self.nucleus_index = int(nucleus_index)
        self.metadata_gate = nn.Sequential(
            nn.Linear(5, max(16, specialist_projection_dim // 4)),
            nn.GELU(),
            nn.Linear(max(16, specialist_projection_dim // 4), 1),
        )
        self.fusion_logit = nn.Parameter(torch.tensor(-1.0))

    def _aggregate_specialist(
        self,
        window_logits: torch.Tensor,
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
        window_valid: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        protein_logits = []
        window_gate = window_logits.new_zeros(window_logits.shape)
        metadata_score = self.metadata_gate(metadata.float()).squeeze(-1)
        for protein_index in range(protein_count):
            selected = (owner == protein_index) & window_valid
            if not selected.any():
                protein_logits.append(window_logits.new_zeros(()))
                continue
            values = window_logits[selected]
            if self.aggregation == "mean":
                aggregate = values.mean()
                window_gate[selected] = 1.0 / values.shape[0]
            elif self.aggregation == "attention":
                gate = torch.softmax(values.float() + metadata_score[selected], dim=0)
                aggregate = torch.sum(gate.to(values.dtype) * values)
                window_gate[selected] = gate.to(window_gate.dtype)
            else:
                probability = torch.sigmoid(values)
                top_count = min(2, probability.shape[0])
                top = torch.topk(probability, top_count, dim=0)
                top_probability = top.values
                selected_gate = torch.zeros_like(probability)
                selected_gate.scatter_(0, top.indices, 1.0 / top_count)
                window_gate[selected] = selected_gate
                noisy_or = 1.0 - torch.prod(1.0 - top_probability)
                aggregate = torch.logit(noisy_or.clamp(1e-6, 1.0 - 1e-6))
            protein_logits.append(aggregate)
        return torch.stack(protein_logits), window_gate

    def forward(
        self,
        windows: torch.Tensor,
        window_mask: torch.Tensor,
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
        window_valid: torch.Tensor | None = None,
    ) -> ModelOutput:
        if window_valid is None:
            window_valid = torch.ones(
                windows.shape[0], dtype=torch.bool, device=windows.device
            )
        general = self.general_mil(
            windows,
            window_mask,
            owner,
            metadata,
            protein_count,
            window_valid=window_valid,
        )
        specialist_window_logit, specialist_attention = self.specialist(
            windows, window_mask
        )
        specialist_logit, specialist_window_gate = self._aggregate_specialist(
            specialist_window_logit,
            owner,
            metadata,
            protein_count,
            window_valid,
        )
        generic_nucleus_logit = general.logits[:, self.nucleus_index]
        fusion_weight = torch.sigmoid(self.fusion_logit)
        logits = general.logits.clone()
        logits[:, self.nucleus_index] = (
            generic_nucleus_logit + fusion_weight * specialist_logit
        )
        auxiliary = dict(general.auxiliary or {})
        auxiliary.update(
            {
                "generic_nucleus_logit": generic_nucleus_logit,
                "specialist_nucleus_logit": specialist_logit,
                "specialist_window_logit": specialist_window_logit,
                "specialist_window_gate": specialist_window_gate,
                "specialist_attention": specialist_attention,
                "fusion_weight": fusion_weight,
            }
        )
        return ModelOutput(
            logits=logits,
            attention=general.attention,
            query_diversity_loss=(
                general.query_diversity_loss
                + query_diversity_loss(self.specialist.query)
            ),
            spectral_loss=(
                general.spectral_loss
                + dct_high_frequency_energy(specialist_attention)
            ),
            evidence=general.evidence,
            auxiliary=auxiliary,
        )


def window_mil_architecture(model_config: dict[str, Any]) -> str:
    architecture = (
        f"{model_config['query_mode']}-{model_config['mil_aggregation']}"
    )
    if bool(model_config.get("nucleus_specialist", False)):
        architecture += "-nucleus-specialist"
    return architecture


def build_window_mil_head(
    input_dim: int,
    output_dim: int,
    nucleus_index: int,
    model_config: dict[str, Any],
) -> nn.Module:
    residue_head = LabelQueryHead(
        input_dim=input_dim,
        projection_dim=int(model_config["projection_dim"]),
        output_dim=output_dim,
        heads=int(model_config.get("heads", 4)),
        mode=model_config["query_mode"],
        dropout=float(model_config.get("dropout", 0.15)),
        smoothing_radius=int(model_config.get("smoothing_radius", 0)),
        dct_cutoff_fraction=float(model_config.get("dct_cutoff_fraction", 0.25)),
    )
    if not bool(model_config.get("nucleus_specialist", False)):
        return WindowMILHead(residue_head, model_config["mil_aggregation"])
    return WindowMILWithNucleusSpecialist(
        residue_head,
        model_config["mil_aggregation"],
        nucleus_index,
        specialist_projection_dim=int(
            model_config.get("specialist_projection_dim", 384)
        ),
        specialist_heads=int(model_config.get("specialist_heads", 4)),
        dropout=float(model_config.get("dropout", 0.15)),
    )


class LabelTokenFusion(nn.Module):
    """Fuse per-label expert tokens without residue-by-residue cross-attention."""

    def __init__(
        self,
        expert_dims: list[int],
        output_dim: int,
        fusion_dim: int = 256,
        heads: int = 4,
        mode: str = "label_token_cross_attention",
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if not 2 <= len(expert_dims) <= 5:
            raise ValueError("Label-token fusion requires 2..5 experts")
        if fusion_dim % heads:
            raise ValueError("fusion_dim must be divisible by heads")
        if mode not in {
            "logit_late_fusion",
            "label_wise_gate",
            "label_token_cross_attention",
            "label_wise_moe",
        }:
            raise ValueError(f"Unsupported fusion mode: {mode}")
        self.output_dim = output_dim
        self.expert_count = len(expert_dims)
        self.fusion_dim = fusion_dim
        self.mode = mode
        self.projections = nn.ModuleList(
            [nn.Sequential(nn.LayerNorm(dim), nn.Linear(dim, fusion_dim)) for dim in expert_dims]
        )
        self.expert_classifiers = nn.Parameter(
            torch.empty(self.expert_count, output_dim, fusion_dim)
        )
        self.expert_bias = nn.Parameter(torch.zeros(self.expert_count, output_dim))
        nn.init.xavier_uniform_(self.expert_classifiers)
        self.static_mixture_logits = nn.Parameter(torch.zeros(output_dim, self.expert_count))
        self.gate = nn.Sequential(
            nn.LayerNorm(self.expert_count * fusion_dim),
            nn.Linear(self.expert_count * fusion_dim, fusion_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, self.expert_count),
        )
        self.label_query = nn.Parameter(torch.empty(output_dim, 1, fusion_dim))
        nn.init.normal_(self.label_query, std=1.0 / math.sqrt(fusion_dim))
        self.cross_attention = nn.MultiheadAttention(
            fusion_dim, heads, dropout=dropout, batch_first=True
        )
        self.fused_classifier = nn.Sequential(
            nn.LayerNorm(fusion_dim),
            nn.Linear(fusion_dim, fusion_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim // 2, 1),
        )

    def forward(
        self,
        expert_tokens: list[torch.Tensor],
        expert_availability: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        if len(expert_tokens) != self.expert_count:
            raise ValueError("Unexpected number of expert token tensors")
        reference_shape = expert_tokens[0].shape[:2]
        if reference_shape[1] != self.output_dim:
            raise ValueError("Expert token label count does not match output_dim")
        projected = []
        for projection, token in zip(self.projections, expert_tokens):
            if token.ndim != 3 or token.shape[:2] != reference_shape:
                raise ValueError("Expert tokens must align as [batch, label, hidden]")
            projected.append(projection(token.float()))
        stacked = torch.stack(projected, dim=2)  # [B, C, E, H]
        batch = stacked.shape[0]
        if expert_availability is None:
            expert_availability = torch.ones(
                (batch, self.expert_count), dtype=torch.bool, device=stacked.device
            )
        else:
            expert_availability = expert_availability.to(
                device=stacked.device, dtype=torch.bool
            )
        if expert_availability.shape != (batch, self.expert_count):
            raise ValueError("expert_availability must have shape [batch, expert]")
        if not expert_availability.any(dim=-1).all():
            raise ValueError("Every protein must have at least one available expert")
        label_availability = expert_availability[:, None, :].expand(
            -1, self.output_dim, -1
        )

        def masked_mixture(logits: torch.Tensor) -> torch.Tensor:
            return torch.softmax(logits.masked_fill(~label_availability, -torch.inf), dim=-1)

        expert_logits = torch.einsum(
            "bceh,ech->bce", stacked, self.expert_classifiers
        ) + self.expert_bias.transpose(0, 1)[None, :, :]

        if self.mode == "logit_late_fusion":
            static_logits = self.static_mixture_logits[None, :, :].expand(batch, -1, -1)
            mixture = masked_mixture(static_logits)
            logits = (expert_logits * mixture).sum(dim=-1)
            return logits, {
                "expert_logits": expert_logits,
                "mixture": mixture,
                "expert_availability": expert_availability,
            }

        dynamic = masked_mixture(self.gate(stacked.flatten(2)))
        if self.mode in {"label_wise_gate", "label_wise_moe"}:
            if self.mode == "label_wise_moe":
                static_logits = self.static_mixture_logits[None, :, :].expand(batch, -1, -1)
                static = masked_mixture(static_logits)
                mixture = 0.5 * dynamic + 0.5 * static
            else:
                mixture = dynamic
            logits = (expert_logits * mixture).sum(dim=-1)
            return logits, {
                "expert_logits": expert_logits,
                "mixture": mixture,
                "expert_availability": expert_availability,
            }

        batch, labels, experts, hidden = stacked.shape
        key_value = stacked.reshape(batch * labels, experts, hidden)
        query = self.label_query[None, :, :, :].expand(batch, -1, -1, -1)
        query = query.reshape(batch * labels, 1, hidden)
        fused, attention = self.cross_attention(
            query,
            key_value,
            key_value,
            key_padding_mask=(~label_availability).reshape(batch * labels, experts),
            need_weights=True,
            average_attn_weights=True,
        )
        fused = fused.reshape(batch, labels, hidden)
        logits = self.fused_classifier(fused).squeeze(-1)
        return logits, {
            "expert_logits": expert_logits,
            "mixture": attention.reshape(batch, labels, experts),
            "fused_tokens": fused,
            "expert_availability": expert_availability,
        }


class MultiBackboneLabelTokenModel(nn.Module):
    """Independent residue pooling per backbone followed by label-token fusion."""

    def __init__(
        self,
        input_dims: list[int],
        projection_dims: list[int],
        output_dim: int,
        fusion_mode: str,
        fusion_dim: int = 256,
        residue_heads: int = 4,
        fusion_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if len(input_dims) != len(projection_dims):
            raise ValueError("input_dims and projection_dims must have equal length")
        self.experts = nn.ModuleList(
            [
                LabelQueryHead(
                    input_dim=input_dim,
                    projection_dim=projection_dim,
                    output_dim=output_dim,
                    heads=residue_heads,
                    mode="label",
                    dropout=dropout,
                )
                for input_dim, projection_dim in zip(input_dims, projection_dims)
            ]
        )
        hidden_dims = [expert.classifier_weight.shape[1] for expert in self.experts]
        self.fusion = LabelTokenFusion(
            hidden_dims,
            output_dim,
            fusion_dim=fusion_dim,
            heads=fusion_heads,
            mode=fusion_mode,
            dropout=dropout,
        )

    def forward(
        self, residues: list[torch.Tensor], masks: list[torch.Tensor]
    ) -> ModelOutput:
        if len(residues) != len(self.experts) or len(masks) != len(self.experts):
            raise ValueError("One residue and mask tensor is required per expert")
        outputs = [
            expert(expert_residues, expert_mask)
            for expert, expert_residues, expert_mask in zip(self.experts, residues, masks)
        ]
        if any(output.evidence is None for output in outputs):
            raise RuntimeError("All fusion experts must emit per-label evidence tokens")
        logits, fusion_auxiliary = self.fusion(
            [output.evidence for output in outputs if output.evidence is not None]
        )
        diversity = torch.stack([output.query_diversity_loss for output in outputs]).mean()
        spectral = torch.stack([output.spectral_loss for output in outputs]).mean()
        minimum_length = min(output.attention.shape[-1] for output in outputs)
        attention = torch.stack(
            [output.attention[..., :minimum_length] for output in outputs], dim=2
        ).mean(dim=2)
        return ModelOutput(
            logits=logits,
            attention=attention,
            query_diversity_loss=diversity,
            spectral_loss=spectral,
            evidence=fusion_auxiliary.get("fused_tokens"),
            auxiliary=fusion_auxiliary,
        )


class MultiBackboneWindowFusion(nn.Module):
    """Window MIL per backbone followed by protein-level label-token fusion."""

    def __init__(
        self,
        input_dims: list[int],
        projection_dims: list[int],
        output_dim: int,
        mil_aggregation: str,
        fusion_mode: str,
        fusion_dim: int = 256,
        heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if len(input_dims) != len(projection_dims):
            raise ValueError("input_dims and projection_dims must have equal length")
        self.experts = nn.ModuleList()
        hidden_dims = []
        for input_dim, projection_dim in zip(input_dims, projection_dims):
            residue_head = LabelQueryHead(
                input_dim,
                projection_dim,
                output_dim,
                heads=heads,
                mode="label",
                dropout=dropout,
            )
            hidden_dims.append(residue_head.classifier_weight.shape[1])
            self.experts.append(WindowMILHead(residue_head, mil_aggregation))
        self.fusion = LabelTokenFusion(
            hidden_dims,
            output_dim,
            fusion_dim=fusion_dim,
            heads=heads,
            mode=fusion_mode,
            dropout=dropout,
        )

    def forward(
        self,
        windows: list[torch.Tensor],
        masks: list[torch.Tensor],
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
        availability: torch.Tensor | None = None,
    ) -> ModelOutput:
        if len(windows) != len(self.experts) or len(masks) != len(self.experts):
            raise ValueError("One window and mask tensor is required per backbone")
        if availability is not None and availability.shape != (
            windows[0].shape[0], windows[0].shape[1], len(self.experts)
        ):
            raise ValueError("availability must have shape [window, residue, expert]")
        outputs = []
        for expert_index, (expert, item, mask) in enumerate(
            zip(self.experts, windows, masks)
        ):
            if availability is None:
                effective_mask = mask
                window_valid = torch.ones(
                    item.shape[0], dtype=torch.bool, device=item.device
                )
                safe_item = item
                safe_mask = mask
            else:
                effective_mask = mask & availability[..., expert_index].bool()
                window_valid = effective_mask.any(dim=-1)
                safe_item = item
                safe_mask = effective_mask.clone()
                invalid_windows = ~window_valid
                if invalid_windows.any():
                    # Avoid all-masked softmax inside the residue head.  These
                    # neutral outputs are excluded by window_valid below.
                    safe_item = item.clone()
                    safe_item[invalid_windows] = 0.0
                    safe_mask[invalid_windows] = mask[invalid_windows]
            outputs.append(
                expert(
                    safe_item,
                    safe_mask,
                    owner,
                    metadata,
                    protein_count,
                    window_valid=window_valid,
                )
            )
        if any(output.evidence is None for output in outputs):
            raise ValueError("Label-token fusion is incompatible with top2_noisy_or window MIL")
        expert_availability = torch.stack(
            [output.auxiliary["protein_availability"] for output in outputs], dim=-1
        )
        logits, auxiliary = self.fusion(
            [output.evidence for output in outputs if output.evidence is not None],
            expert_availability=expert_availability,
        )
        return ModelOutput(
            logits=logits,
            attention=outputs[0].attention,
            query_diversity_loss=torch.stack(
                [output.query_diversity_loss for output in outputs]
            ).mean(),
            spectral_loss=torch.stack([output.spectral_loss for output in outputs]).mean(),
            evidence=auxiliary.get("fused_tokens"),
            auxiliary=auxiliary,
        )


class ScalarMixedWindowMIL(nn.Module):
    """Learn a normalized scalar mixture of aligned layers before window MIL."""

    def __init__(
        self,
        input_dim: int,
        layer_count: int,
        projection_dim: int,
        output_dim: int,
        mil_aggregation: str,
        heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if layer_count < 2:
            raise ValueError("Scalar mixing requires at least two encoder layers")
        self.layer_norms = nn.ModuleList([nn.LayerNorm(input_dim) for _ in range(layer_count)])
        self.scalar_logits = nn.Parameter(torch.zeros(layer_count))
        self.gamma = nn.Parameter(torch.ones(()))
        residue_head = LabelQueryHead(
            input_dim,
            projection_dim,
            output_dim,
            heads=heads,
            mode="label",
            dropout=dropout,
        )
        self.head = WindowMILHead(residue_head, mil_aggregation)

    def forward(
        self,
        windows: list[torch.Tensor],
        masks: list[torch.Tensor],
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
    ) -> ModelOutput:
        if len(windows) != len(self.layer_norms) or len(masks) != len(self.layer_norms):
            raise ValueError("Scalar mixing needs one aligned tensor per frozen layer")
        reference_shape = windows[0].shape
        reference_mask = masks[0]
        for window, mask in zip(windows[1:], masks[1:]):
            if window.shape != reference_shape or not torch.equal(mask, reference_mask):
                raise ValueError("Scalar-mixed layers must have identical window/residue alignment")
        weights = torch.softmax(self.scalar_logits, dim=0)
        normalized = [normalizer(window.float()) for normalizer, window in zip(self.layer_norms, windows)]
        mixed = self.gamma * torch.stack(normalized, dim=0).mul(
            weights[:, None, None, None]
        ).sum(dim=0)
        output = self.head(mixed, reference_mask, owner, metadata, protein_count)
        auxiliary = dict(output.auxiliary or {})
        auxiliary["layer_mixture"] = weights
        auxiliary["layer_mixture_gamma"] = self.gamma
        output.auxiliary = auxiliary
        return output


class AlignedResidueGatedFusion(nn.Module):
    """Fuse position-aligned AA/3Di residues with an explicit availability mask."""

    def __init__(
        self,
        input_dims: list[int],
        projection_dim: int,
        output_dim: int,
        heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if len(input_dims) < 2:
            raise ValueError("Aligned residue fusion requires at least two backbones")
        self.projections = nn.ModuleList(
            [nn.Linear(input_dim, projection_dim) for input_dim in input_dims]
        )
        self.gate = nn.Sequential(
            nn.LayerNorm(projection_dim * len(input_dims) + len(input_dims)),
            nn.Linear(projection_dim * len(input_dims) + len(input_dims), projection_dim),
            nn.GELU(),
            nn.Linear(projection_dim, len(input_dims)),
        )
        self.head = LabelQueryHead(
            projection_dim, projection_dim, output_dim, heads=heads, mode="label", dropout=dropout
        )

    def forward(
        self,
        residues: list[torch.Tensor],
        masks: list[torch.Tensor],
        availability: torch.Tensor | None = None,
    ) -> ModelOutput:
        if len(residues) != len(self.projections) or len(masks) != len(self.projections):
            raise ValueError("One aligned residue and mask tensor is required per backbone")
        shape = residues[0].shape[:2]
        if any(tensor.shape[:2] != shape for tensor in residues):
            raise ValueError("Aligned residue fusion requires equal batch and residue dimensions")
        projected = torch.stack(
            [projection(tensor.float()) for projection, tensor in zip(self.projections, residues)],
            dim=2,
        )
        mask_stack = torch.stack(masks, dim=-1)
        if availability is None:
            availability = mask_stack
        availability = availability.bool() & mask_stack
        if not torch.all(availability.any(dim=-1)):
            raise ValueError("Every residue must be available from at least one expert")
        gate_input = torch.cat(
            (projected.flatten(2), availability.to(projected.dtype)), dim=-1
        )
        gate_logits = self.gate(gate_input).masked_fill(
            ~availability, torch.finfo(projected.dtype).min
        )
        gate = torch.softmax(gate_logits.float(), dim=-1).to(projected.dtype)
        fused = (projected * gate[..., None]).sum(dim=2)
        output = self.head(fused, availability.any(dim=-1))
        auxiliary = dict(output.auxiliary or {})
        auxiliary["residue_expert_gate"] = gate
        output.auxiliary = auxiliary
        return output


class AlignedWindowResidueFusion(nn.Module):
    """Position-aligned residue gating followed by label-wise window MIL."""

    def __init__(
        self,
        input_dims: list[int],
        projection_dim: int,
        output_dim: int,
        aggregation: str = "attention",
        heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if aggregation not in {"mean", "attention", "top2_noisy_or"}:
            raise ValueError(f"Unsupported MIL aggregation: {aggregation}")
        self.residue_fusion = AlignedResidueGatedFusion(
            input_dims, projection_dim, output_dim, heads=heads, dropout=dropout
        )
        self.aggregation = aggregation
        evidence_dim = int(self.residue_fusion.head.classifier_weight.shape[1])
        self.metadata_projection = nn.Sequential(
            nn.Linear(5, evidence_dim),
            nn.GELU(),
            nn.Linear(evidence_dim, evidence_dim),
        )
        self.gate_weight = nn.Parameter(torch.empty(output_dim, evidence_dim))
        self.gate_bias = nn.Parameter(torch.zeros(output_dim))
        nn.init.xavier_uniform_(self.gate_weight)

    def classify(self, evidence: torch.Tensor) -> torch.Tensor:
        head = self.residue_fusion.head
        return (
            torch.einsum("bch,ch->bc", evidence, head.classifier_weight)
            + head.classifier_bias
        )

    def forward(
        self,
        windows: list[torch.Tensor],
        masks: list[torch.Tensor],
        owner: torch.Tensor,
        metadata: torch.Tensor,
        protein_count: int,
        availability: torch.Tensor,
    ) -> ModelOutput:
        window_output = self.residue_fusion(windows, masks, availability)
        if window_output.evidence is None:
            raise RuntimeError("Aligned window MIL requires per-label evidence tokens")
        evidence = window_output.evidence + self.metadata_projection(metadata.float())[:, None, :]
        protein_evidence = []
        protein_logits = []
        for protein_index in range(protein_count):
            selected_evidence = evidence[owner == protein_index]
            if selected_evidence.shape[0] == 0:
                raise ValueError(f"Protein {protein_index} has no windows")
            if self.aggregation == "mean":
                protein_evidence.append(selected_evidence.mean(dim=0))
            elif self.aggregation == "attention":
                gate_logits = torch.einsum(
                    "wch,ch->wc", selected_evidence, self.gate_weight
                ) + self.gate_bias
                gate = torch.softmax(gate_logits, dim=0)
                protein_evidence.append(
                    torch.einsum("wc,wch->ch", gate, selected_evidence)
                )
            else:
                probability = torch.sigmoid(self.classify(selected_evidence))
                top_probability = torch.topk(
                    probability, min(2, probability.shape[0]), dim=0
                ).values
                noisy_or = 1.0 - torch.prod(1.0 - top_probability, dim=0)
                protein_logits.append(torch.logit(noisy_or.clamp(1e-6, 1.0 - 1e-6)))
        if self.aggregation == "top2_noisy_or":
            logits = torch.stack(protein_logits)
            aggregated_evidence = None
        else:
            aggregated_evidence = torch.stack(protein_evidence)
            logits = self.classify(aggregated_evidence)
        auxiliary = dict(window_output.auxiliary or {})
        auxiliary["structure_availability_fraction"] = availability.float().mean()
        return ModelOutput(
            logits=logits,
            attention=window_output.attention,
            query_diversity_loss=window_output.query_diversity_loss,
            spectral_loss=window_output.spectral_loss,
            evidence=aggregated_evidence,
            auxiliary=auxiliary,
        )


class SparseConceptExpert(nn.Module):
    """Small expert over selected SAE latents and fixed per-protein aggregates."""

    def __init__(
        self,
        feature_count: int,
        aggregate_count: int,
        output_dim: int,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if not 1 <= feature_count <= 128:
            raise ValueError("InterPLM concept expert allows 1..128 selected latents")
        self.feature_count = feature_count
        self.aggregate_count = aggregate_count
        self.output_dim = output_dim
        self.encoder = nn.Sequential(
            nn.LayerNorm(feature_count * aggregate_count),
            nn.Linear(feature_count * aggregate_count, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
        )
        self.label_tokens = nn.Parameter(torch.empty(output_dim, hidden_dim))
        self.label_bias = nn.Parameter(torch.zeros(output_dim))
        nn.init.xavier_uniform_(self.label_tokens)

    def forward(self, aggregates: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if aggregates.shape[1:] != (self.feature_count, self.aggregate_count):
            raise ValueError("Concept aggregates must be [batch, selected_latent, statistic]")
        hidden = self.encoder(aggregates.float().flatten(1))
        logits = hidden @ self.label_tokens.transpose(0, 1) + self.label_bias
        return logits, hidden


class ConceptAttentionPrior(nn.Module):
    """Map selected residue-level SAE activations to additive label attention priors."""

    def __init__(self, feature_count: int, output_dim: int, initial_scale: float = 0.1) -> None:
        super().__init__()
        self.feature_count = feature_count
        self.output_dim = output_dim
        self.weight = nn.Parameter(torch.zeros(output_dim, feature_count))
        self.raw_scale = nn.Parameter(torch.tensor(math.log(math.expm1(initial_scale))))

    def forward(self, activation: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        if activation.ndim != 3 or activation.shape[-1] != self.feature_count:
            raise ValueError("SAE activation must be [batch, residue, selected_latent]")
        prior = torch.einsum("blk,ck->bcl", activation.float(), self.weight)
        prior = prior.masked_fill(~mask[:, None, :], 0.0)
        scale = torch.nn.functional.softplus(self.raw_scale)
        return scale * prior


class GeneralWithConceptAttentionPrior(nn.Module):
    """Inject residue-level sparse concepts into label-query attention logits."""

    def __init__(
        self,
        general_head: LabelQueryHead,
        feature_count: int,
        initial_scale: float = 0.1,
    ) -> None:
        super().__init__()
        self.general_head = general_head
        self.concept_prior = ConceptAttentionPrior(
            feature_count, general_head.output_dim, initial_scale=initial_scale
        )

    def forward(
        self,
        residues: torch.Tensor,
        mask: torch.Tensor,
        concept_activation: torch.Tensor,
    ) -> ModelOutput:
        prior = self.concept_prior(concept_activation, mask)
        output = self.general_head(residues, mask, attention_prior=prior)
        auxiliary = dict(output.auxiliary or {})
        auxiliary.update(
            {
                "concept_attention_prior_scale": torch.nn.functional.softplus(
                    self.concept_prior.raw_scale
                ),
                "concept_attention_prior_mean_abs": prior.abs().mean(),
            }
        )
        return ModelOutput(
            logits=output.logits,
            attention=output.attention,
            query_diversity_loss=output.query_diversity_loss,
            spectral_loss=output.spectral_loss,
            evidence=output.evidence,
            auxiliary=auxiliary,
        )


class GeneralWithConceptExpert(nn.Module):
    """Fuse a general residue model with a sparse InterPLM concept expert."""

    def __init__(
        self,
        general_head: LabelQueryHead,
        concept_expert: SparseConceptExpert,
        mode: str = "independent_expert",
    ) -> None:
        super().__init__()
        if mode not in {"independent_expert", "specialist_token"}:
            raise ValueError(f"Unsupported concept fusion mode: {mode}")
        self.general_head = general_head
        self.concept_expert = concept_expert
        self.mode = mode
        self.label_gate = nn.Parameter(torch.full((general_head.output_dim,), -1.0))
        general_hidden = general_head.classifier_weight.shape[1]
        concept_hidden = concept_expert.label_tokens.shape[1]
        self.concept_to_label = nn.Linear(concept_hidden, general_hidden)

    def forward(
        self,
        residues: torch.Tensor,
        mask: torch.Tensor,
        concept_aggregates: torch.Tensor,
    ) -> ModelOutput:
        general = self.general_head(residues, mask)
        concept_logits, concept_hidden = self.concept_expert(concept_aggregates)
        gate = torch.sigmoid(self.label_gate)[None, :]
        if self.mode == "independent_expert":
            logits = general.logits + gate * concept_logits
            evidence = general.evidence
        else:
            if general.evidence is None:
                raise RuntimeError("specialist_token mode requires general label evidence")
            token = self.concept_to_label(concept_hidden)[:, None, :]
            evidence = general.evidence + gate[:, :, None] * token
            logits = torch.einsum(
                "bch,ch->bc", evidence, self.general_head.classifier_weight
            ) + self.general_head.classifier_bias
        auxiliary = dict(general.auxiliary or {})
        auxiliary.update(
            {"concept_logits": concept_logits, "concept_gate": gate, "concept_hidden": concept_hidden}
        )
        return ModelOutput(
            logits=logits,
            attention=general.attention,
            query_diversity_loss=general.query_diversity_loss,
            spectral_loss=general.spectral_loss,
            evidence=evidence,
            auxiliary=auxiliary,
        )


class AsymmetricLoss(nn.Module):
    def __init__(
        self,
        gamma_negative: float = 4.0,
        gamma_positive: float = 1.0,
        clip: float = 0.05,
        label_smoothing: float = 0.0,
        positive_weight: torch.Tensor | None = None,
    ):
        super().__init__()
        self.gamma_negative = gamma_negative
        self.gamma_positive = gamma_positive
        self.clip = clip
        self.label_smoothing = label_smoothing
        self.register_buffer("positive_weight", positive_weight.float() if positive_weight is not None else None)

    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        if self.label_smoothing:
            target = target * (1.0 - self.label_smoothing) + 0.5 * self.label_smoothing
        probability = torch.sigmoid(logits)
        positive = probability
        negative = 1.0 - probability
        if self.clip:
            negative = (negative + self.clip).clamp(max=1.0)
        loss = target * torch.log(positive.clamp_min(1e-8))
        loss += (1.0 - target) * torch.log(negative.clamp_min(1e-8))
        correct_probability = positive * target + negative * (1.0 - target)
        weight = torch.pow(
            (1.0 - correct_probability).clamp_min(0.0),
            self.gamma_positive * target + self.gamma_negative * (1.0 - target),
        )
        loss = -(loss * weight)
        if self.positive_weight is not None:
            loss = loss * self.positive_weight.view(1, -1)
        return loss.mean()


class WeightedFocalLoss(nn.Module):
    def __init__(self, positive_weight: torch.Tensor, gamma: float = 2.0):
        super().__init__()
        self.register_buffer("positive_weight", positive_weight.float())
        self.gamma = gamma

    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        base = F.binary_cross_entropy_with_logits(logits, target, reduction="none")
        probability = torch.sigmoid(logits)
        correct_probability = probability * target + (1.0 - probability) * (1.0 - target)
        class_weight = target * self.positive_weight + (1.0 - target)
        return (base * torch.pow(1.0 - correct_probability, self.gamma) * class_weight).mean()


class DistributionBalancedFocalLoss(nn.Module):
    """Distribution-balanced focal loss with repeat-rate and negative tolerance."""

    def __init__(
        self,
        label_frequency: torch.Tensor,
        sample_count: int,
        focal_gamma: float = 2.0,
        rebalance_beta: float = 10.0,
        rebalance_alpha: float = 0.1,
        rebalance_gamma: float = 0.2,
        negative_scale: float = 5.0,
        init_bias_factor: float = 0.05,
    ):
        super().__init__()
        frequency = label_frequency.float().clamp_min(1.0)
        inverse = 1.0 / frequency
        prior = (frequency / sample_count).clamp(1e-6, 1.0 - 1e-6)
        initial_bias = -torch.log((1.0 - prior) / prior) * init_bias_factor
        self.register_buffer("inverse_frequency", inverse)
        self.register_buffer("initial_bias", initial_bias)
        self.focal_gamma = focal_gamma
        self.rebalance_beta = rebalance_beta
        self.rebalance_alpha = rebalance_alpha
        self.rebalance_gamma = rebalance_gamma
        self.negative_scale = negative_scale

    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        repeat_rate = (target * self.inverse_frequency).sum(dim=1, keepdim=True).clamp_min(1e-8)
        positive_repeat_weight = self.inverse_frequency[None, :] / repeat_rate
        rebalanced = (
            torch.sigmoid(
                self.rebalance_beta * (positive_repeat_weight - self.rebalance_alpha)
            )
            + self.rebalance_gamma
        )
        weight = 1.0 + target * (rebalanced - 1.0)
        regulated = logits - self.initial_bias[None, :]
        regulated = regulated * (target + (1.0 - target) * self.negative_scale)
        probability = torch.sigmoid(regulated)
        base = F.binary_cross_entropy_with_logits(regulated, target, reduction="none")
        focal = torch.pow(torch.abs(target - probability), self.focal_gamma)
        return (base * focal * weight).mean()


def build_loss(
    name: str,
    label_frequency: torch.Tensor,
    sample_count: int,
    options: dict[str, Any] | None = None,
) -> nn.Module:
    options = options or {}
    if name == "asl":
        positive_weight = None
        if options.get("use_positive_weight", False):
            negative = sample_count - label_frequency.float()
            positive_weight = (negative / label_frequency.float().clamp_min(1.0)).clamp(
                max=float(options.get("max_positive_weight", 20.0))
            )
        return AsymmetricLoss(
            gamma_negative=float(options.get("gamma_negative", 4.0)),
            gamma_positive=float(options.get("gamma_positive", 1.0)),
            clip=float(options.get("clip", 0.05)),
            label_smoothing=float(options.get("label_smoothing", 0.0)),
            positive_weight=positive_weight,
        )
    if name == "deeploc_focal":
        if label_frequency.numel() != 11:
            raise ValueError("DeepLoc fixed focal weights are defined only for the 11-label task")
        return WeightedFocalLoss(
            torch.tensor([1, 1, 1, 3, 2.3, 4, 9.5, 4.5, 6.6, 7.7, 32], dtype=torch.float32)
        )
    if name == "sqrt_focal":
        positive_weight = torch.sqrt(sample_count / label_frequency.float().clamp_min(1.0))
        return WeightedFocalLoss(positive_weight / positive_weight.mean())
    if name == "distribution_balanced_focal":
        return DistributionBalancedFocalLoss(label_frequency, sample_count)
    raise ValueError(f"Unknown loss: {name}")


def sorting_signal_kl_loss(
    attention: torch.Tensor,
    supervision: torch.Tensor,
    supervised_labels: torch.Tensor,
) -> torch.Tensor:
    """KL(target||attention) for allowed training-fold sorting annotations only.

    attention: [B, C, L], supervision: [B, C, L], supervised_labels: [B, C].
    """
    if attention.shape != supervision.shape:
        raise ValueError("Sorting-signal attention and supervision shapes differ")
    if supervised_labels.shape != attention.shape[:-1]:
        raise ValueError("Sorting-signal label mask shape differs from attention")

    # Attention is emitted in the autocast dtype.  With fp16, 1e-8 rounds to
    # zero, so log(attention.clamp_min(1e-8)) can still be -inf.  Multiplying
    # that value by a zero target (especially for unsupervised specialist
    # queries) produces NaN and corrupts the first optimizer step.  Evaluate
    # the regularizer in fp32 and explicitly exclude unsupervised channels.
    attention_fp32 = attention.float()
    supervision_fp32 = supervision.float()
    weights = supervised_labels.float()
    active = weights > 0
    if not active.any():
        return attention_fp32.sum() * 0.0
    if not torch.isfinite(attention_fp32).all():
        raise FloatingPointError("Sorting-signal attention contains nonfinite values")
    if not torch.isfinite(supervision_fp32).all():
        raise FloatingPointError("Sorting-signal supervision contains nonfinite values")
    if (supervision_fp32[active].sum(dim=-1) <= 0).any():
        raise ValueError("A supervised sorting-signal channel has zero target mass")

    target = supervision_fp32 / supervision_fp32.sum(
        dim=-1, keepdim=True
    ).clamp_min(1e-8)
    log_target = torch.log(target.clamp_min(1e-8))
    log_attention = torch.log(attention_fp32.clamp_min(1e-8))
    per_label = (target * (log_target - log_attention)).sum(dim=-1)
    per_label = torch.where(active, per_label, torch.zeros_like(per_label))
    return (per_label * weights).sum() / weights.sum().clamp_min(1.0)


def conflict_repulsion_loss(
    attention: torch.Tensor,
    conflict_mask: torch.Tensor,
    negative_labels: torch.Tensor,
) -> torch.Tensor:
    """Discourage attention on conflicting signal residues only for negative labels."""
    mass = (attention * conflict_mask).sum(dim=-1)
    weights = negative_labels.float()
    return (mass * weights).sum() / weights.sum().clamp_min(1.0)
