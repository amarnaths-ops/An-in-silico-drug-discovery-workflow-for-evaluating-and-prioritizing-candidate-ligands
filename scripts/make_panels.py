import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =====================================================
# SETTINGS - change things here only
# =====================================================
RANKING_FILE = "results/ranking.csv"
CLUSTER_FILE = "data/haddock_clusters.csv"
IMAGE_DIR = "images"
OUT_DIR = "results/panels"
ONLY_TOP = True      # True = only top compounds, False = all compounds
TOP_N = 5
# Image names expected in the images folder (png, jpg or jpeg):
#   <compound_id>_haddock.png   original HADDOCK plot
#   <compound_id>_pymol.png     PyMOL picture
#   <compound_id>_2d.png        Protein Plus 2D interactions
# =====================================================

os.makedirs(OUT_DIR, exist_ok=True)
rank = pd.read_csv(RANKING_FILE).sort_values("weighted_rank")
if ONLY_TOP:
    rank = rank.head(TOP_N)

cl = pd.read_csv(CLUSTER_FILE) if os.path.exists(CLUSTER_FILE) else None
if cl is not None:
    cl["compound_id"] = cl["compound_id"].astype(str)


def find_image(compound, suffix):
    for ext in (".png", ".jpg", ".jpeg"):
        path = os.path.join(IMAGE_DIR, compound + suffix + ext)
        if os.path.exists(path):
            return path
    return None


def show_image(ax, compound, suffix, title):
    ax.set_title(title, fontsize=11)
    path = find_image(compound, suffix)
    if path is None:
        ax.text(0.5, 0.5, "Image not found:\n" + compound + suffix + ".png",
                ha="center", va="center", fontsize=10, color="red")
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.imshow(plt.imread(path))
        ax.axis("off")


def model_plot(ax, compound):
    ax.set_title("Model generated: HADDOCK score vs iRMSD", fontsize=11)
    if cl is None or compound not in set(cl["compound_id"]):
        ax.text(0.5, 0.5, "No cluster data", ha="center", va="center")
        return
    # other compounds faded in the background for context
    others = cl[cl["compound_id"] != compound]
    ax.scatter(others["irmsd"], others["haddock_score"], s=15,
               color="lightgrey", label="other compounds")
    sub = cl[cl["compound_id"] == compound]
    ax.errorbar(sub["irmsd"], sub["haddock_score"],
                xerr=sub["irmsd_sd"].fillna(0), yerr=sub["haddock_score_sd"].fillna(0),
                fmt="none", ecolor="steelblue", elinewidth=1)
    ax.scatter(sub["irmsd"], sub["haddock_score"], s=sub["cluster_size"] * 3,
               color="steelblue", edgecolors="black", label=compound)
    ax.set_xlabel("iRMSD (A)")
    ax.set_ylabel("HADDOCK score")
    ax.legend(fontsize=7, loc="best")


for _, row in rank.iterrows():
    cid = str(row["compound_id"])
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.5))
    model_plot(axes[0], cid)
    show_image(axes[1], cid, "_haddock", "Original HADDOCK plot")
    show_image(axes[2], cid, "_pymol", "PyMOL")
    show_image(axes[3], cid, "_2d", "Protein Plus 2D interactions")
    fig.suptitle(cid + "   (rank " + str(int(row["weighted_rank"])) +
                 ", score " + format(row["weighted_score"], ".2f") + ")",
                 fontsize=14)
    plt.tight_layout()
    out = os.path.join(OUT_DIR, cid + "_panel.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved:", out)
