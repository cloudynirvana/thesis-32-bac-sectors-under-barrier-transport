#!/usr/bin/env python3
"""Lumped versus barrier-coupled BAC sector maps on one five-index toy.

The weight matrix, the grounded Laplacian, and the sector classifier are the
objects of Thesis #25 (itself the classifier of Thesis #18). Each off-diagonal
pair, under the barrier reading, is a desmoplastic edge in the sense of
Thesis #21: a shedding rate and a conductance, collapsed by the Schur
complement to the series flux φ = λ κ / (λ + κ).

The lumped map decays the weight itself. The barrier map decays shedding and
holds the conductance schedule. A co-decay check multiplies both rates by the
same exponential, which restores the lumped weight.

Seed 20260921 is used only for the monotonicity draws. The sector grids are
deterministic. Regenerating this file rewrites sim/results.json and
sim/figures/.

Research only. Not a medical device. Not a dose. Not a rejuvenation claim.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921

LABELS = ["molecular", "cellular", "tissue", "organism", "evolutionary"]
N = 5
ORG = 3
GROUND = 4

DELTA_AGE = 0.04
DELTA_CUT = 0.08
DELTA_AMB_CUT = 0.06
DELTA_AMB_BLOCK = 0.02
SIGMA_HOLD = 0.15
T_END = 80.0
N_FINE = 401

THR_LAM_FALL = 0.5
THR_BOTH_FALL = 0.5
THR_SECTOR_BAND = 1.5
THR_CUT_COLLAPSE = 0.25
THR_BLOCK_HELD = 0.9
THR_SECTOR_RISE = 3.0

N_GRID = 13
# Stalled fractions of Thesis #21, recomputed from that deposit's (λ, κ).
LAM_T21 = np.array([0.09, 0.06, 0.04])
KAP_T21 = np.array([0.25, 0.025, 0.04])
S_T21 = LAM_T21 / (LAM_T21 + KAP_T21)
S_SHED = float(S_T21[0])
S_TRANSPORT = float(S_T21[1])
S_BALANCED = float(S_T21[2])
S_TRANSPARENT = 1.0e-12

KINDS = ("held", "aging", "cancer", "ambiguous", "diagonal")


def initial_weights() -> np.ndarray:
    pairs = {
        (0, 1): 0.80,
        (0, 2): 0.35,
        (0, 3): 0.55,
        (0, 4): 0.20,
        (1, 2): 0.75,
        (1, 3): 0.60,
        (1, 4): 0.25,
        (2, 3): 0.70,
        (2, 4): 0.30,
        (3, 4): 0.40,
    }
    w = np.zeros((N, N), dtype=float)
    for (i, j), value in pairs.items():
        w[i, j] = w[j, i] = value
    return w


def edge_index() -> list[tuple[int, int, str]]:
    edges = []
    for i in range(N):
        for j in range(i + 1, N):
            sector = "cut" if (i == ORG or j == ORG) else "block"
            edges.append((i, j, sector))
    return edges


EDGES = edge_index()


def symmetrize(weights: np.ndarray) -> np.ndarray:
    a = 0.5 * (np.asarray(weights, dtype=float) + np.asarray(weights, dtype=float).T)
    np.fill_diagonal(a, 0.0)
    return a


def laplacian(weights: np.ndarray) -> np.ndarray:
    a = symmetrize(weights)
    return np.diag(a.sum(axis=1)) - a


def algebraic_connectivity(weights: np.ndarray) -> float:
    ev = np.linalg.eigvalsh(laplacian(weights))
    return float(np.sort(ev)[1])


def grounded_lambda_min(weights: np.ndarray, ground: int = GROUND) -> float:
    ell = laplacian(weights)
    keep = [i for i in range(N) if i != ground]
    return float(np.linalg.eigvalsh(ell[np.ix_(keep, keep)])[0])


def cut_and_block_means(weights: np.ndarray) -> tuple[float, float]:
    a = symmetrize(weights)
    cut = []
    block = []
    for i, j, sector in EDGES:
        if sector == "cut":
            cut.append(a[i, j])
        else:
            block.append(a[i, j])
    return float(np.mean(cut)), float(np.mean(block))


def edge_rates(kind: str) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for i, j, sector in EDGES:
        if kind == "held":
            rate = 0.0
        elif kind == "aging":
            rate = DELTA_AGE
        elif kind == "cancer":
            rate = DELTA_CUT if sector == "cut" else 0.0
        elif kind == "ambiguous":
            rate = DELTA_AMB_CUT if sector == "cut" else DELTA_AMB_BLOCK
        elif kind == "diagonal":
            rate = 0.0
        else:
            raise ValueError(kind)
        rates[i, j] = rates[j, i] = rate
    return rates


def plane_rates(delta_cut: float, delta_block: float) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for i, j, sector in EDGES:
        rate = delta_cut if sector == "cut" else delta_block
        rates[i, j] = rates[j, i] = rate
    return rates


def series_factor(delta: float, s: float, t: float) -> float:
    """φ(t) / W(0) when λ(t) = λ(0) e^{−δ t} and κ is held, with s = λ(0)/(λ(0)+κ)."""
    if delta == 0.0 or t == 0.0:
        return 1.0
    e = float(np.exp(-delta * t))
    return e / (s * e + (1.0 - s))


def weight_matrix(t: float, rates: np.ndarray, schedule: dict, w0: np.ndarray) -> np.ndarray:
    w = np.zeros((N, N), dtype=float)
    kind = schedule["kind"]
    for i, j, sector in EDGES:
        delta = float(rates[i, j])
        base = float(w0[i, j])
        if kind == "lumped":
            value = base * float(np.exp(-delta * t))
        elif kind == "codecay":
            # Both rates carry δ. The stalled fraction is held at 1/2 only as a
            # split; Proposition 3 says the split cancels.
            s = 0.5
            e = float(np.exp(-delta * t))
            lam = (base / (1.0 - s)) * e
            kap = (base / s) * e
            value = explicit_phi(lam, kap)
        elif kind == "barrier":
            s = float(schedule["s_cut"] if sector == "cut" else schedule["s_block"])
            value = base * series_factor(delta, s, t)
        else:
            raise ValueError(kind)
        w[i, j] = w[j, i] = value
    return w


def classify(lam_ratio: float, block_ratio: float, cut_ratio: float, sector_change: float) -> str:
    aging = (
        lam_ratio < THR_LAM_FALL
        and block_ratio < THR_BOTH_FALL
        and cut_ratio < THR_BOTH_FALL
        and (1.0 / THR_SECTOR_BAND) < sector_change < THR_SECTOR_BAND
    )
    cancer = (
        lam_ratio < THR_LAM_FALL
        and cut_ratio < THR_CUT_COLLAPSE
        and block_ratio > THR_BLOCK_HELD
        and sector_change > THR_SECTOR_RISE
    )
    if aging and cancer:
        return "both"
    if aging:
        return "aging-like"
    if cancer:
        return "cancer-like"
    return "neither"


def sector_record(rates: np.ndarray, schedule: dict, w0: np.ndarray, t_end: float = T_END) -> dict:
    w_a = weight_matrix(0.0, rates, schedule, w0)
    w_b = weight_matrix(t_end, rates, schedule, w0)
    lam0 = grounded_lambda_min(w_a)
    lam1 = grounded_lambda_min(w_b)
    cut0, block0 = cut_and_block_means(w_a)
    cut1, block1 = cut_and_block_means(w_b)
    lam_ratio = lam1 / lam0
    cut_ratio = cut1 / cut0
    block_ratio = block1 / block0
    sector0 = block0 / cut0
    sector1 = block1 / cut1
    sector_change = sector1 / sector0
    return {
        "class": classify(lam_ratio, block_ratio, cut_ratio, sector_change),
        "lambda_min_0": lam0,
        "lambda_min_T": lam1,
        "lambda_min_ratio": lam_ratio,
        "mean_cut_0": cut0,
        "mean_cut_T": cut1,
        "mean_cut_ratio": cut_ratio,
        "mean_block_0": block0,
        "mean_block_T": block1,
        "mean_block_ratio": block_ratio,
        "sector_change": sector_change,
        "lambda2_0": algebraic_connectivity(w_a),
        "lambda2_T": algebraic_connectivity(w_b),
    }


def lambda_path(rates: np.ndarray, schedule: dict, w0: np.ndarray, t: np.ndarray) -> np.ndarray:
    out = np.empty(t.size, dtype=float)
    for k, time in enumerate(t):
        out[k] = grounded_lambda_min(weight_matrix(float(time), rates, schedule, w0))
    return out


def first_crossing(t: np.ndarray, lam: np.ndarray, sigma: float) -> float | None:
    v = lam - sigma
    if np.all(v > 0):
        return None
    idx = int(np.argmax(v <= 0))
    if idx == 0:
        return float(t[0])
    t0, t1 = float(t[idx - 1]), float(t[idx])
    v0, v1 = float(v[idx - 1]), float(v[idx])
    if v1 == v0:
        return t1
    return t0 + (0.0 - v0) * (t1 - t0) / (v1 - v0)


def schedules() -> list[dict]:
    return [
        {"name": "lumped", "kind": "lumped", "s_cut": None, "s_block": None},
        {"name": "codecay", "kind": "codecay", "s_cut": None, "s_block": None},
        {
            "name": "transparent",
            "kind": "barrier",
            "s_cut": S_TRANSPARENT,
            "s_block": S_TRANSPARENT,
        },
        {
            "name": "uniform_shedding",
            "kind": "barrier",
            "s_cut": S_SHED,
            "s_block": S_SHED,
        },
        {
            "name": "uniform_balanced",
            "kind": "barrier",
            "s_cut": S_BALANCED,
            "s_block": S_BALANCED,
        },
        {
            "name": "uniform_transport",
            "kind": "barrier",
            "s_cut": S_TRANSPORT,
            "s_block": S_TRANSPORT,
        },
        {
            "name": "cut_transport",
            "kind": "barrier",
            "s_cut": S_TRANSPORT,
            "s_block": S_SHED,
        },
        {
            "name": "block_transport",
            "kind": "barrier",
            "s_cut": S_SHED,
            "s_block": S_TRANSPORT,
        },
    ]


def set_compare(left: list[str], right: list[str], label: str) -> dict:
    a = np.array([x == label for x in left], dtype=bool)
    b = np.array([x == label for x in right], dtype=bool)
    both = int(np.sum(a & b))
    left_only = int(np.sum(a & ~b))
    right_only = int(np.sum(~a & b))
    union = both + left_only + right_only
    return {
        "label": label,
        "n": int(a.size),
        "n_lumped": int(np.sum(a)),
        "n_barrier": int(np.sum(b)),
        "both": both,
        "lumped_only": left_only,
        "barrier_only": right_only,
        "sets_equal": bool(left_only == 0 and right_only == 0),
        "jaccard": None if union == 0 else both / union,
    }


def label_agreement(left: list[str], right: list[str]) -> dict:
    same = int(sum(a == b for a, b in zip(left, right)))
    n = len(left)
    return {"n": n, "n_agree": same, "n_differ": n - same, "fraction_agree": same / n}


def counts(labels: list[str]) -> dict:
    return {
        "aging-like": int(sum(x == "aging-like" for x in labels)),
        "cancer-like": int(sum(x == "cancer-like" for x in labels)),
        "neither": int(sum(x == "neither" for x in labels)),
        "both": int(sum(x == "both" for x in labels)),
    }


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )


PATH_COLORS = {
    "held": "#4C5760",
    "aging": "#1F4E79",
    "cancer": "#8C3A3A",
    "ambiguous": "#8A6A12",
    "diagonal": "#2F6B4F",
}
PATH_NAMES = {
    "held": "held weights",
    "aging": "global decay",
    "cancer": "organism-scale cut",
    "ambiguous": "uneven decay",
    "diagonal": "diagonal only",
}
CLASS_COLOR = {
    "aging-like": "#1F4E79",
    "cancer-like": "#8C3A3A",
    "neither": "#C8C2B4",
    "both": "#6B4C9A",
}


def _scatter_classes(ax, xs, ys, labels, title: str) -> None:
    for lab, color in CLASS_COLOR.items():
        pts = [(x, y) for x, y, z in zip(xs, ys, labels) if z == lab]
        if not pts:
            continue
        ax.scatter(
            [p[0] for p in pts],
            [p[1] for p in pts],
            s=36,
            c=color,
            linewidths=0.0,
            zorder=2,
        )
    ax.set_title(title)
    ax.set_xlabel("block rate")
    ax.set_xlim(-0.004, DELTA_CUT + 0.004)
    ax.set_ylim(-0.004, DELTA_CUT + 0.004)


def draw(paths: dict, planes: dict, s_plane: dict, deltas: np.ndarray) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    t = paths["t"]

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.9), sharey=True)
    for ax, name, title in (
        (axes[0], "lumped", "lumped coupling"),
        (axes[1], "cut_transport", "cut transport-limited"),
    ):
        for key in KINDS:
            ax.plot(t, paths[name][key], color=PATH_COLORS[key], lw=1.6, label=PATH_NAMES[key])
        ax.axhline(SIGMA_HOLD, color="#666666", ls=":", lw=0.8)
        ax.set_xlabel("toy time")
        ax.set_title(title)
    axes[0].set_ylabel("λ_min of the grounded operator")
    axes[1].legend(frameon=False, fontsize=7.2, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "lambda_min_lumped_vs_barrier.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.45), sharex=True, sharey=True)
    panel = (
        ("lumped", "lumped coupling"),
        ("cut_transport", "cut transport-limited"),
        ("block_transport", "block transport-limited"),
    )
    for ax, (name, title) in zip(axes, panel):
        block = planes[name]
        _scatter_classes(
            ax,
            [p["delta_block"] for p in block],
            [p["delta_cut"] for p in block],
            [p["class"] for p in block],
            title,
        )
        ax.plot([0.0, DELTA_CUT], [0.0, DELTA_CUT], color="#888888", lw=0.6, zorder=1)
    axes[0].set_ylabel("cut rate")
    fig.legend(
        handles=[
            Patch(facecolor=CLASS_COLOR["aging-like"], edgecolor="none", label="aging-like"),
            Patch(facecolor=CLASS_COLOR["cancer-like"], edgecolor="none", label="cancer-like"),
            Patch(facecolor=CLASS_COLOR["neither"], edgecolor="none", label="neither"),
        ],
        frameon=False,
        fontsize=8,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(FIG / "rate_plane_classes.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), sharex=True, sharey=True)
    for ax, name, title in (
        (axes[0], "cut_transport", "disagreement, cut transport-limited"),
        (axes[1], "block_transport", "disagreement, block transport-limited"),
    ):
        lumped = { (p["delta_cut"], p["delta_block"]): p["class"] for p in planes["lumped"] }
        for p in planes[name]:
            key = (p["delta_cut"], p["delta_block"])
            same = p["class"] == lumped[key]
            ax.scatter(
                p["delta_block"],
                p["delta_cut"],
                s=36,
                c="#1F4E79" if same else "#8C3A3A",
                linewidths=0.0,
            )
        ax.plot([0.0, DELTA_CUT], [0.0, DELTA_CUT], color="#888888", lw=0.6, zorder=0)
        ax.set_title(title)
        ax.set_xlabel("block rate")
        ax.set_xlim(-0.004, DELTA_CUT + 0.004)
        ax.set_ylim(-0.004, DELTA_CUT + 0.004)
    axes[0].set_ylabel("cut rate")
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#1F4E79", markersize=6, label="same class"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#8C3A3A", markersize=6, label="different class"),
        ],
        frameon=False,
        fontsize=8,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(FIG / "rate_plane_disagreement.png")
    plt.close(fig)

    s_values = np.array(s_plane["s_values"], dtype=float)
    labels = s_plane["global_decay"]
    code = {"neither": 0.0, "aging-like": 1.0, "cancer-like": 2.0, "both": 3.0}
    # Stored rows vary the block stalled fraction; columns vary the cut fraction.
    # imshow uses the first axis as y, so transpose puts the cut fraction on y.
    grid = np.array([[code[lab] for lab in row] for row in labels], dtype=float).T
    cmap = matplotlib.colors.ListedColormap(
        [CLASS_COLOR["neither"], CLASS_COLOR["aging-like"], CLASS_COLOR["cancer-like"], CLASS_COLOR["both"]]
    )
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ax.imshow(
        grid,
        origin="lower",
        cmap=cmap,
        vmin=-0.5,
        vmax=3.5,
        extent=(
            s_values[0] - 0.5 * (s_values[1] - s_values[0]),
            s_values[-1] + 0.5 * (s_values[1] - s_values[0]),
            s_values[0] - 0.5 * (s_values[1] - s_values[0]),
            s_values[-1] + 0.5 * (s_values[1] - s_values[0]),
        ),
        interpolation="nearest",
        aspect="equal",
    )
    ax.plot([float(s_values[0]), float(s_values[-1])], [float(s_values[0]), float(s_values[-1])], color="#222222", lw=0.7)
    ax.scatter([S_SHED], [S_TRANSPORT], s=46, c="white", edgecolors="#111111", linewidths=0.8, zorder=3)
    ax.scatter([S_TRANSPORT], [S_SHED], s=46, marker="s", c="white", edgecolors="#111111", linewidths=0.8, zorder=3)
    ax.set_xlabel("block stalled fraction")
    ax.set_ylabel("cut stalled fraction")
    ax.set_title("global decay, barrier schedule")
    fig.legend(
        handles=[
            Patch(facecolor=CLASS_COLOR["aging-like"], edgecolor="none", label="aging-like"),
            Patch(facecolor=CLASS_COLOR["cancer-like"], edgecolor="none", label="cancer-like"),
            Patch(facecolor=CLASS_COLOR["neither"], edgecolor="none", label="neither"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#111111", markersize=6, label="cut transport-limited"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor="white", markeredgecolor="#111111", markersize=6, label="block transport-limited"),
        ],
        frameon=False,
        fontsize=7.4,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    fig.savefig(FIG / "barrier_schedule_plane.png")
    plt.close(fig)


def round_floats(obj, nd=12):
    if isinstance(obj, float):
        if obj != obj or obj in (float("inf"), float("-inf")):
            return None
        return round(obj, nd)
    if isinstance(obj, dict):
        return {k: round_floats(v, nd) for k, v in obj.items()}
    if isinstance(obj, list):
        return [round_floats(v, nd) for v in obj]
    if isinstance(obj, (np.floating,)):
        return round_floats(float(obj), nd)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def explicit_phi(lam: float, kap: float) -> float:
    return lam * kap / (lam + kap)


def main() -> None:
    w0 = initial_weights()
    scheds = {s["name"]: s for s in schedules()}
    t_fine = np.linspace(0.0, T_END, N_FINE)
    deltas = np.linspace(0.0, DELTA_CUT, N_GRID)
    s_values = np.linspace(0.05, 0.95, N_GRID)

    named = {}
    for name, schedule in scheds.items():
        named[name] = {}
        for kind in KINDS:
            named[name][kind] = sector_record(edge_rates(kind), schedule, w0)

    plane_names = ("lumped", "uniform_transport", "cut_transport", "block_transport", "transparent", "codecay")
    planes = {name: [] for name in plane_names}
    for delta_block in deltas:
        for delta_cut in deltas:
            rates = plane_rates(float(delta_cut), float(delta_block))
            for name in plane_names:
                rec = sector_record(rates, scheds[name], w0)
                planes[name].append(
                    {
                        "delta_cut": float(delta_cut),
                        "delta_block": float(delta_block),
                        "class": rec["class"],
                        "lambda_min_ratio": rec["lambda_min_ratio"],
                        "mean_cut_ratio": rec["mean_cut_ratio"],
                        "mean_block_ratio": rec["mean_block_ratio"],
                        "sector_change": rec["sector_change"],
                    }
                )

    s_plane = {"s_values": s_values.tolist(), "global_decay": [], "organism_cut": [], "uneven": []}
    rate_for_s = {
        "global_decay": edge_rates("aging"),
        "organism_cut": edge_rates("cancer"),
        "uneven": edge_rates("ambiguous"),
    }
    for key, rates in rate_for_s.items():
        rows = []
        for s_block in s_values:
            row = []
            for s_cut in s_values:
                schedule = {
                    "name": "grid",
                    "kind": "barrier",
                    "s_cut": float(s_cut),
                    "s_block": float(s_block),
                }
                row.append(sector_record(rates, schedule, w0)["class"])
            rows.append(row)
        s_plane[key] = rows

    paths = {"t": t_fine}
    for name in ("lumped", "cut_transport", "block_transport", "uniform_transport"):
        paths[name] = {
            kind: lambda_path(edge_rates(kind), scheds[name], w0, t_fine) for kind in KINDS
        }
    crossings = {}
    for name in ("lumped", "cut_transport"):
        crossings[name] = {
            kind: first_crossing(t_fine, paths[name][kind], SIGMA_HOLD) for kind in KINDS
        }

    # --- identities, fixed before the maps are treated as a claim ---
    rng = np.random.default_rng(SEED)
    mono_fail = 0
    for _ in range(200):
        i, j = rng.choice(N, size=2, replace=False)
        lifted = w0.copy()
        lifted[i, j] += 0.05
        lifted[j, i] += 0.05
        if grounded_lambda_min(lifted) + 1e-10 < grounded_lambda_min(w0):
            mono_fail += 1

    ray_err = 0.0
    for _ in range(50):
        x = rng.normal(size=N)
        x = x - x.mean()
        x = x / np.linalg.norm(x)
        quad = float(x @ laplacian(w0) @ x)
        acc = 0.0
        a = symmetrize(w0)
        for i in range(N):
            for j in range(i + 1, N):
                acc += a[i, j] * (x[i] - x[j]) ** 2
        ray_err = max(ray_err, abs(quad - acc))

    phi_err = 0.0
    for s in (S_SHED, S_BALANCED, S_TRANSPORT, 0.1, 0.9):
        for delta in (0.0, DELTA_AGE, DELTA_CUT):
            for t in (0.0, 20.0, 80.0):
                w = 0.55
                lam0 = w / (1.0 - s)
                kap = w / s
                lam = lam0 * float(np.exp(-delta * t))
                direct = explicit_phi(lam, kap)
                via = w * series_factor(delta, s, t)
                phi_err = max(phi_err, abs(direct - via))

    schur_err = 0.0
    for lam, kap in ((0.8, 0.3), (0.2, 0.9), (0.06, 0.025)):
        phi = explicit_phi(lam, kap)
        leaves = np.array([[lam, 0.0], [0.0, kap]])
        link = np.array([-lam, -kap])
        schur = leaves - np.outer(link, link) / (lam + kap)
        expect = phi * np.array([[1.0, -1.0], [-1.0, 1.0]])
        schur_err = max(schur_err, float(np.max(np.abs(schur - expect))))

    expected = {
        "held": "neither",
        "aging": "aging-like",
        "cancer": "cancer-like",
        "ambiguous": "neither",
        "diagonal": "neither",
    }
    lumped_ok = {k: named["lumped"][k]["class"] == expected[k] for k in KINDS}
    codecay_match = all(named["codecay"][k]["class"] == named["lumped"][k]["class"] for k in KINDS)
    transparent_match = all(
        named["transparent"][k]["class"] == named["lumped"][k]["class"] for k in KINDS
    )
    codecay_gap = 0.0
    transparent_gap = 0.0
    for kind in KINDS:
        for key in ("lambda_min_ratio", "mean_cut_ratio", "mean_block_ratio", "sector_change"):
            codecay_gap = max(codecay_gap, abs(named["codecay"][kind][key] - named["lumped"][kind][key]))
            transparent_gap = max(
                transparent_gap, abs(named["transparent"][kind][key] - named["lumped"][kind][key])
            )

    uniform_sector_err = 0.0
    for uname in ("uniform_shedding", "uniform_balanced", "uniform_transport", "transparent"):
        uniform_sector_err = max(uniform_sector_err, abs(named[uname]["aging"]["sector_change"] - 1.0))

    aging_l = named["lumped"]["aging"]
    scale = aging_l["lambda_min_0"] * float(np.exp(-DELTA_AGE * T_END))
    aging_scale_err = abs(aging_l["lambda_min_T"] - scale)
    analytic_cross = float(np.log(aging_l["lambda_min_0"] / SIGMA_HOLD) / DELTA_AGE)
    cross_err = abs(float(crossings["lumped"]["aging"]) - analytic_cross)

    cut0, block0 = cut_and_block_means(w0)
    both_count = 0
    for name, rows in planes.items():
        both_count += sum(p["class"] == "both" for p in rows)
    for name in scheds:
        both_count += sum(named[name][k]["class"] == "both" for k in KINDS)
    for key in ("global_decay", "organism_cut", "uneven"):
        both_count += sum(lab == "both" for row in s_plane[key] for lab in row)

    plane_compare = {}
    lumped_labels = [p["class"] for p in planes["lumped"]]
    for name in plane_names:
        if name == "lumped":
            continue
        labels = [p["class"] for p in planes[name]]
        plane_compare[name] = {
            "agreement": label_agreement(lumped_labels, labels),
            "aging": set_compare(lumped_labels, labels, "aging-like"),
            "cancer": set_compare(lumped_labels, labels, "cancer-like"),
            "counts_lumped": counts(lumped_labels),
            "counts_barrier": counts(labels),
        }

    named_compare = {}
    for name in scheds:
        if name == "lumped":
            continue
        left = [named["lumped"][k]["class"] for k in KINDS]
        right = [named[name][k]["class"] for k in KINDS]
        named_compare[name] = {
            "agreement": label_agreement(left, right),
            "by_kind": {k: {"lumped": named["lumped"][k]["class"], "barrier": named[name][k]["class"]} for k in KINDS},
        }

    s_global_flat = [lab for row in s_plane["global_decay"] for lab in row]
    s_summary = {
        "global_decay": counts(s_global_flat),
        "organism_cut": counts([lab for row in s_plane["organism_cut"] for lab in row]),
        "uneven": counts([lab for row in s_plane["uneven"] for lab in row]),
    }
    # Diagonal of the stalled-fraction plane is the uniform schedule, at global decay.
    diag_classes = [s_plane["global_decay"][i][i] for i in range(N_GRID)]

    checks = {
        "monotonicity_grounded_failures": mono_fail,
        "rayleigh_max_abs_err": float(ray_err),
        "series_flux_max_abs_err": float(phi_err),
        "schur_edge_max_abs_err": float(schur_err),
        "lumped_classes_match_predeclaration": lumped_ok,
        "codecay_classes_match_lumped": codecay_match,
        "transparent_classes_match_lumped": transparent_match,
        "codecay_ratio_gap": float(codecay_gap),
        "transparent_ratio_gap": float(transparent_gap),
        "uniform_aging_sector_change_err": float(uniform_sector_err),
        "aging_scale_max_abs_err": float(aging_scale_err),
        "aging_analytic_crossing": analytic_cross,
        "aging_numeric_crossing_abs_err": float(cross_err),
        "initial_lambda_min": float(named["lumped"]["held"]["lambda_min_0"]),
        "initial_lambda2": float(named["lumped"]["held"]["lambda2_0"]),
        "initial_mean_cut": cut0,
        "initial_mean_block": block0,
        "both_class_count": both_count,
        "n_edges": len(EDGES),
        "n_cut": sum(1 for *_, s in EDGES if s == "cut"),
        "n_block": sum(1 for *_, s in EDGES if s == "block"),
    }

    hard = [
        checks["monotonicity_grounded_failures"] == 0,
        checks["rayleigh_max_abs_err"] < 1e-9,
        checks["series_flux_max_abs_err"] < 1e-12,
        checks["schur_edge_max_abs_err"] < 1e-12,
        all(lumped_ok.values()),
        checks["codecay_classes_match_lumped"] is True,
        checks["transparent_classes_match_lumped"] is True,
        checks["codecay_ratio_gap"] < 1e-12,
        checks["transparent_ratio_gap"] < 1e-9,
        checks["uniform_aging_sector_change_err"] < 1e-12,
        checks["aging_scale_max_abs_err"] < 1e-10,
        checks["aging_numeric_crossing_abs_err"] < 1e-3,
        abs(checks["initial_lambda_min"] - 0.285102) < 5e-6,
        abs(checks["initial_lambda2"] - 1.413084) < 5e-6,
        abs(checks["initial_mean_cut"] - 0.5625) < 1e-12,
        abs(checks["initial_mean_block"] - 0.4416666667) < 1e-9,
        checks["both_class_count"] == 0,
        checks["n_cut"] == 4 and checks["n_block"] == 6,
        # Heterogeneous T21 split moves the global-decay sector ratio off 1.
        abs(named["cut_transport"]["aging"]["sector_change"] - 1.0) > 0.2,
        abs(named["block_transport"]["aging"]["sector_change"] - 1.0) > 0.2,
        # And off the aging band, while the lumped path stays inside it.
        named["lumped"]["aging"]["class"] == "aging-like",
        named["cut_transport"]["aging"]["class"] == "neither",
        named["block_transport"]["aging"]["class"] == "neither",
        not ((1.0 / THR_SECTOR_BAND) < named["cut_transport"]["aging"]["sector_change"] < THR_SECTOR_BAND),
        not ((1.0 / THR_SECTOR_BAND) < named["block_transport"]["aging"]["sector_change"] < THR_SECTOR_BAND),
        plane_compare["codecay"]["agreement"]["n_differ"] == 0,
        plane_compare["transparent"]["agreement"]["n_differ"] == 0,
        plane_compare["cut_transport"]["aging"]["sets_equal"] is False,
        plane_compare["block_transport"]["aging"]["sets_equal"] is False,
    ]
    if not all(hard):
        raise SystemExit("toy checks failed:\n" + json.dumps(round_floats(checks), indent=2))

    # Strip time series from the JSON; figures already hold the curves.
    named_out = {}
    for name, by_kind in named.items():
        named_out[name] = {}
        for kind, rec in by_kind.items():
            row = dict(rec)
            row["crossing_hold"] = crossings.get(name, {}).get(kind)
            named_out[name][kind] = row

    diagonal = []
    by_key = {
        name: {(p["delta_cut"], p["delta_block"]): p for p in planes[name]}
        for name in ("lumped", "cut_transport", "block_transport", "uniform_transport")
    }
    for delta in deltas:
        key = (float(delta), float(delta))
        diagonal.append(
            {
                "delta": float(delta),
                "lumped": by_key["lumped"][key]["class"],
                "uniform_transport": by_key["uniform_transport"][key]["class"],
                "cut_transport": by_key["cut_transport"][key]["class"],
                "block_transport": by_key["block_transport"][key]["class"],
                "lumped_sector_change": by_key["lumped"][key]["sector_change"],
                "cut_sector_change": by_key["cut_transport"][key]["sector_change"],
                "block_sector_change": by_key["block_transport"][key]["sector_change"],
                "cut_lambda_ratio": by_key["cut_transport"][key]["lambda_min_ratio"],
                "lumped_lambda_ratio": by_key["lumped"][key]["lambda_min_ratio"],
            }
        )

    payload = {
        "seed": SEED,
        "labels": LABELS,
        "organism_index": ORG,
        "ground_index": GROUND,
        "depends_on": ["T25", "T21"],
        "parents_of_parents": ["T18", "T13", "T11", "T05"],
        "shared_horizon": {"t0": 0.0, "t_end": T_END, "n_fine": N_FINE},
        "stalled_fractions": {
            "shedding_limited": S_SHED,
            "transport_limited": S_TRANSPORT,
            "balanced": S_BALANCED,
            "source_lambda": LAM_T21.tolist(),
            "source_kappa": KAP_T21.tolist(),
            "transparent": S_TRANSPARENT,
        },
        "parameters": {
            "delta_age": DELTA_AGE,
            "delta_cut": DELTA_CUT,
            "delta_ambiguous_cut": DELTA_AMB_CUT,
            "delta_ambiguous_block": DELTA_AMB_BLOCK,
            "sigma_hold": SIGMA_HOLD,
            "n_grid": N_GRID,
        },
        "thresholds": {
            "lambda_ratio_fall": THR_LAM_FALL,
            "both_sectors_fall": THR_BOTH_FALL,
            "sector_change_band": THR_SECTOR_BAND,
            "cut_collapse": THR_CUT_COLLAPSE,
            "block_held": THR_BLOCK_HELD,
            "sector_rise": THR_SECTOR_RISE,
        },
        "schedules": schedules(),
        "initial": {
            "lambda_min": checks["initial_lambda_min"],
            "lambda2": checks["initial_lambda2"],
            "mean_cut": cut0,
            "mean_block": block0,
        },
        "named": named_out,
        "named_compare": named_compare,
        "plane": {
            "deltas": deltas.tolist(),
            "compare": plane_compare,
            "counts": {name: counts([p["class"] for p in rows]) for name, rows in planes.items()},
            "equal_rate_diagonal": diagonal,
        },
        "barrier_plane": {
            "s_values": s_values.tolist(),
            "counts": s_summary,
            "global_decay_diagonal_classes": diag_classes,
            "classes": {
                "global_decay": s_plane["global_decay"],
                "organism_cut": s_plane["organism_cut"],
                "uneven": s_plane["uneven"],
            },
        },
        "checks": checks,
        "witnesses": {
            "mild_cancer_nodes": [],
        },
    }
    for dc in (5.0 / 150.0, 6.0 / 150.0):
        def nearest(name: str, cut: float, block: float) -> dict:
            best = min(
                planes[name],
                key=lambda p: abs(p["delta_cut"] - cut) + abs(p["delta_block"] - block),
            )
            return best
        lum = nearest("lumped", dc, 0.0)
        cutb = nearest("cut_transport", dc, 0.0)
        payload["witnesses"]["mild_cancer_nodes"].append(
            {
                "delta_cut": lum["delta_cut"],
                "delta_block": lum["delta_block"],
                "lumped_class": lum["class"],
                "lumped_lambda_ratio": lum["lambda_min_ratio"],
                "cut_transport_class": cutb["class"],
                "cut_transport_lambda_ratio": cutb["lambda_min_ratio"],
            }
        )

    draw(
        {"t": t_fine, "lumped": paths["lumped"], "cut_transport": paths["cut_transport"]},
        planes,
        {"s_values": s_values.tolist(), "global_decay": s_plane["global_decay"]},
        deltas,
    )

    out = ROOT / "results.json"
    out.write_text(json.dumps(round_floats(payload), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print(
        f"lambda_min0={checks['initial_lambda_min']:.6f} "
        f"lambda2={checks['initial_lambda2']:.6f} "
        f"cross={analytic_cross:.4f}"
    )
    print("--- named classes ---")
    header = f"{'schedule':22} " + " ".join(f"{k[:6]:>8}" for k in KINDS)
    print(header)
    for name in scheds:
        cols = " ".join(f"{named[name][k]['class'][:8]:>8}" for k in KINDS)
        print(f"{name:22} {cols}")
    print("--- global decay ratios (r_λ, r_cut, r_block, s) ---")
    for name in scheds:
        r = named[name]["aging"]
        print(
            f"{name:22} {r['lambda_min_ratio']:.6f} {r['mean_cut_ratio']:.6f} "
            f"{r['mean_block_ratio']:.6f} {r['sector_change']:.6f} {r['class']}"
        )
    print("--- organism-cut ratios ---")
    for name in scheds:
        r = named[name]["cancer"]
        print(
            f"{name:22} {r['lambda_min_ratio']:.6f} {r['mean_cut_ratio']:.6f} "
            f"{r['mean_block_ratio']:.6f} {r['sector_change']:.6f} {r['class']}"
        )
    print("--- uneven ratios ---")
    for name in scheds:
        r = named[name]["ambiguous"]
        print(
            f"{name:22} {r['lambda_min_ratio']:.6f} {r['mean_cut_ratio']:.6f} "
            f"{r['mean_block_ratio']:.6f} {r['sector_change']:.6f} {r['class']}"
        )
    print("--- plane ---")
    for name, comp in plane_compare.items():
        ag = comp["aging"]
        ca = comp["cancer"]
        print(
            f"{name:22} differ={comp['agreement']['n_differ']:3} "
            f"aging L/B/both/Lonly/Bonly={ag['n_lumped']}/{ag['n_barrier']}/{ag['both']}/{ag['lumped_only']}/{ag['barrier_only']} "
            f"cancer L/B/both={ca['n_lumped']}/{ca['n_barrier']}/{ca['both']} "
            f"jA={ag['jaccard']} jC={ca['jaccard']}"
        )
    print("--- s plane counts ---")
    for key, val in s_summary.items():
        print(key, val)
    print("diagonal", diag_classes)
    print("crossings lumped", crossings["lumped"])
    print("crossings cut", crossings["cut_transport"])


if __name__ == "__main__":
    main()
