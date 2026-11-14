import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # saves plots without needing a display
import matplotlib.pyplot as plt

# =====================================================
# SETTINGS - change things here only
# =====================================================
DATA_FILE = "data/main_data.csv"
CLUSTER_FILE = "data/haddock_clusters.csv"
TOP_N = 5        # use 5 for your real 12 compounds
RUNS = 1000      # number of random weight sets for the stability test

directions = {
    "vina_kcal": "lower",
    "cbdock_score": "lower",
    "haddock_score": "lower",
    "prodigy_dg_kcal": "lower",
    "plip_saltbridge_no": "higher",
    "plip_hydrogen_no": "higher",
    "plip_hydrophobic_no": "higher",
    "vina_contacts_4A": "higher",
    "vina_min_dist_pocket": "lower",
}

groups = {
    "docking": ["vina_kcal", "cbdock_score"],
    "haddock": ["haddock_score"],
    "affinity": ["prodigy_dg_kcal"],
    "interactions": ["plip_saltbridge_no", "plip_hydrogen_no", "plip_hydrophobic_no"],
    "site": ["vina_contacts_4A", "vina_min_dist_pocket"],
}

weights = {"docking": 0.25, "haddock": 0.25, "affinity": 0.15,
           "interactions": 0.15, "site": 0.20}
# =====================================================

os.makedirs("results", exist_ok=True)
df = pd.read_csv(DATA_FILE)

# ---------- 1. CHECKS ----------
problems = []

if "compound_id" not in df.columns:
    print("STOP. Column 'compound_id' is missing in", DATA_FILE)
    raise SystemExit
df["compound_id"] = df["compound_id"].astype(str)

dups = df[df["compound_id"].duplicated()]["compound_id"].tolist()
if dups:
    problems.append("Duplicate compound IDs: " + str(dups))

for col in directions:
    if col not in df.columns:
        problems.append("Column missing in CSV: " + col)
        continue

    # keep only the number before "+/-" (or "±"), e.g. "-7.5 +/- 0.3" -> -7.5
    text = (df[col].astype(str)
            .str.replace("±", "+/-", regex=False)
            .str.split("+/-", regex=False).str[0]
            .str.strip())
    numbers = pd.to_numeric(text.where(df[col].notna()), errors="coerce")

    bad = df[numbers.isna() & df[col].notna()]["compound_id"].tolist()
    if bad:
        problems.append("Not a number in " + col + " for: " + str(bad))
    df[col] = numbers

    if numbers.nunique() < 2:
        problems.append("Column has no variation (or is empty): " + col)

    if directions[col] == "lower":
        neg = (numbers < 0).sum()
        pos = (numbers > 0).sum()
        if neg >= 0.75 * numbers.notna().sum() and pos > 0:
            odd = df[numbers > 0]["compound_id"].tolist()
            print("WARNING:", col, "is mostly negative but positive for", odd)

for g, cols in groups.items():
    for c in cols:
        if c not in directions:
            problems.append("Column " + c + " is in groups but not in directions")
if set(groups) != set(weights):
    problems.append("groups and weights do not have the same names")

if problems:
    print("STOP. Fix these problems first:")
    for p in problems:
        print(" -", p)
    raise SystemExit

# ---------- 2. SCALING (0 to 1, 1 = best) ----------
for col, d in directions.items():
    low = df[col].min()
    high = df[col].max()
    if d == "lower":
        df[col + "_norm"] = (high - df[col]) / (high - low)
    else:
        df[col + "_norm"] = (df[col] - low) / (high - low)

for g, cols in groups.items():
    df[g] = df[[c + "_norm" for c in cols]].mean(axis=1)


def consensus(use):
    """Weighted mean of group scores; missing groups are skipped and
    the remaining weights are re-scaled."""
    total = 0
    wsum = 0
    for g in use:
        total = total + df[g].fillna(0) * weights[g]
        wsum = wsum + df[g].notna() * weights[g]
    return total / wsum.replace(0, np.nan)


# ---------- 3. RANKING (two methods) ----------
df["weighted_score"] = consensus(groups)
df["weighted_rank"] = df["weighted_score"].rank(
    ascending=False, method="min", na_option="bottom").astype(int)
df["mean_rank"] = df[list(groups)].rank(ascending=False).mean(axis=1)
df["mean_rank_order"] = df["mean_rank"].rank(
    method="min", na_option="bottom").astype(int)
df["groups_missing"] = len(groups) - df[list(groups)].notna().sum(axis=1)

# ---------- 4. STABILITY TEST ----------
n_top = min(TOP_N, len(df))
rng = np.random.default_rng(42)
names = list(groups)
base = np.array([weights[g] for g in names])
scores = df[names].to_numpy(dtype=float)
present = ~np.isnan(scores)
counts = np.zeros(len(df))

for _ in range(RUNS):
    w = rng.dirichlet(base * 20)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (np.nan_to_num(scores) * w).sum(axis=1) / (present * w).sum(axis=1)
    c = np.nan_to_num(c, nan=-1.0)  # compounds with no data go to the bottom
    counts[np.argsort(-c)[:n_top]] += 1

df["stability_pct"] = (100 * counts / RUNS).round(0)

# ---------- 5. PRINT RESULTS ----------
df = df.sort_values(["weighted_rank", "mean_rank_order"]).reset_index(drop=True)
cols = ["compound_id", "weighted_score", "weighted_rank", "mean_rank_order",
        "stability_pct", "groups_missing"]
print("RANKING (all compounds):")
print(df[cols].round(2).to_string(index=False))

print()
print("TOP", n_top, ":", df["compound_id"].head(n_top).tolist())

print()
print("Leave one group out (rank of each compound):")
for g in groups:
    use = [x for x in groups if x != g]
    r = consensus(use).rank(ascending=False, method="min",
                            na_option="bottom").astype(int)
    print("without", g, ":", dict(zip(df["compound_id"], r)))

print()
for _, row in df.iterrows():
    if row["groups_missing"] > 0:
        print("WARNING:", row["compound_id"], "is missing", int(row["groups_missing"]),
              "group(s). Check before trusting its rank.")
    if row["weighted_rank"] != row["mean_rank_order"]:
        print("NOTE:", row["compound_id"], "has different ranks in the two methods.")

df[cols].round(3).to_csv("results/ranking.csv", index=False)

# ---------- 6. PLOTS ----------
fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(df["compound_id"].tolist()[::-1],
        df["weighted_score"].tolist()[::-1], color="steelblue")
ax.set_xlabel("Consensus score (1 = best)")
ax.set_title("Consensus ranking")
plt.tight_layout()
plt.savefig("results/bar_chart.png", dpi=200)
plt.close()

data = df[names].to_numpy(dtype=float)
fig, ax = plt.subplots(figsize=(7, 4))
im = ax.imshow(data, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=30, ha="right")
ax.set_yticks(range(len(df)))
ax.set_yticklabels(df["compound_id"])
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        if not np.isnan(data[i, j]):
            ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center", color="black")
fig.colorbar(im, label="Score (1 = best)")
ax.set_title("Group scores per compound")
plt.tight_layout()
plt.savefig("results/heatmap.png", dpi=200)
plt.close()

needed = ["compound_id", "irmsd", "irmsd_sd", "haddock_score",
          "haddock_score_sd", "cluster_size"]
missing_cols = []
SHOW_Z = False   # set True to print Z-scores next to top compounds

if os.path.exists(CLUSTER_FILE):
    cl = pd.read_csv(CLUSTER_FILE)
    missing_cols = [c for c in needed if c not in cl.columns]
    if missing_cols:
        print("Cluster file is missing columns", missing_cols,
              "- skipping the HADDOCK vs iRMSD plot.")
    else:
        from matplotlib.lines import Line2D
        cl["compound_id"] = cl["compound_id"].astype(str)
        has_z = SHOW_Z and "z_score" in cl.columns
        top = df["compound_id"].head(n_top).tolist()
        colors = list(plt.cm.tab10.colors) + list(plt.cm.Dark2.colors)
        fig, ax = plt.subplots(figsize=(9, 6))
        handles = []
        for k, (cid, sub) in enumerate(cl.groupby("compound_id")):
            is_top = cid in top
            alpha = 0.9 if is_top else 0.25
            col = colors[k % len(colors)]
            ax.errorbar(sub["irmsd"], sub["haddock_score"],
                        xerr=sub["irmsd_sd"].fillna(0),
                        yerr=sub["haddock_score_sd"].fillna(0),
                        fmt="none", ecolor=col, alpha=alpha * 0.6, elinewidth=0.8)
            ax.scatter(sub["irmsd"], sub["haddock_score"],
                       s=sub["cluster_size"] * 3, color=col, alpha=alpha,
                       edgecolors="black", linewidths=0.6)
            if has_z and is_top:
                for _, r in sub.iterrows():
                    if pd.notna(r["z_score"]):
                        ax.annotate("Z=" + format(r["z_score"], ".1f"),
                                    (r["irmsd"], r["haddock_score"]),
                                    textcoords="offset points", xytext=(6, 6),
                                    fontsize=7)
            handles.append(Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor=col, markeredgecolor="black",
                                  markersize=7, alpha=1.0 if is_top else 0.4,
                                  label=cid + (" (top)" if is_top else "")))
        ax.set_xlabel("iRMSD from best-scoring model (A)")
        ax.set_ylabel("HADDOCK score")
        ax.set_title("HADDOCK score vs iRMSD (dot size = cluster size)")
        ax.legend(handles=handles, fontsize=7, loc="center left",
                  bbox_to_anchor=(1.01, 0.5))
        plt.tight_layout()
        plt.savefig("results/haddock_vs_irmsd.png", dpi=200, bbox_inches="tight")
        plt.close()
else:
    print("No cluster file found, skipping the HADDOCK vs iRMSD plot.")

print()
saved = "ranking.csv, bar_chart.png, heatmap.png"
if os.path.exists(CLUSTER_FILE) and not missing_cols:
    saved += ", haddock_vs_irmsd.png"
print("Saved in results/:", saved)
