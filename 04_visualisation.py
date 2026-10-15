"""Create the final MIC-BPC project figures using Matplotlib."""

import matplotlib

# Use a non-interactive backend so the figures can be saved from the terminal.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Read the cleaned dataset and the results produced in R.

data = pd.read_csv("pseudomonas_MIC_BPC_with_ratio.csv")

results = pd.read_csv("correlation_results.csv").iloc[0]

bootstrap_values = pd.read_csv(
    "bootstrap_correlation_values.csv"
)["spearman_rho"]


# 2. Check that the required data are available.

required_columns = {
    "isolate_id",
    "MIC_tobramycin_ug_ml",
    "BPC_tobramycin_ug_ml",
}

if not required_columns.issubset(data.columns):
    raise ValueError(
        "The cleaned dataset does not contain the required columns."
    )

if len(data) != 30:
    raise ValueError(
        "The cleaned dataset should contain 30 isolates."
    )

if bootstrap_values.empty:
    raise ValueError(
        "The bootstrap correlation file does not contain any values."
    )


# 3. Prepare the data for Figure 1.

# Count how many isolates share each MIC-BPC combination.

point_counts = (
    data.groupby(
        [
            "MIC_tobramycin_ug_ml",
            "BPC_tobramycin_ug_ml",
        ]
    )
    .size()
    .reset_index(name="number_of_isolates")
)


# The concentrations follow doubling steps.
# Log2 spacing places these concentration levels at equal distances.

point_counts["MIC_log2"] = np.log2(
    point_counts["MIC_tobramycin_ug_ml"]
)

point_counts["BPC_log2"] = np.log2(
    point_counts["BPC_tobramycin_ug_ml"]
)


# 4. Create Figure 1: association between MIC and BPC.

# The second panel is reserved for the statistical summary.
# This prevents the summary from covering the data points.

fig = plt.figure(figsize=(12, 7))

layout = fig.add_gridspec(
    1,
    2,
    width_ratios=[4, 1.5],
    wspace=0.15,
)

ax = fig.add_subplot(layout[0, 0])

information_panel = fig.add_subplot(layout[0, 1])

information_panel.axis("off")


# Create the bubble plot.

scatter = ax.scatter(
    point_counts["MIC_log2"],
    point_counts["BPC_log2"],
    s=80 + point_counts["number_of_isolates"] * 75,
    c=point_counts["number_of_isolates"],
    cmap="viridis",
    edgecolor="black",
    alpha=0.85,
)


# Write the number of isolates inside bubbles representing repeated values.

for _, row in point_counts.iterrows():

    if row["number_of_isolates"] > 1:

        ax.text(
            row["MIC_log2"],
            row["BPC_log2"],
            str(int(row["number_of_isolates"])),
            horizontalalignment="center",
            verticalalignment="center",
            color="white",
            fontweight="bold",
        )


# Obtain the original MIC and BPC concentration levels.

mic_ticks = np.sort(
    data["MIC_tobramycin_ug_ml"].unique()
)

bpc_ticks = np.sort(
    data["BPC_tobramycin_ug_ml"].unique()
)


# Display the original concentrations as axis labels.

ax.set_xticks(
    np.log2(mic_ticks),
    labels=[str(value) for value in mic_ticks],
)

ax.set_yticks(
    np.log2(bpc_ticks),
    labels=[str(value) for value in bpc_ticks],
)

# Add the title and axis labels.
ax.set_title(
    "Association Between Tobramycin MIC and BPC",
    fontweight="bold",
    pad=15
)

ax.set_ylabel(
    "Tobramycin BPC (µg/mL)",
    labelpad=10
)
ax.set_xlabel(
    "Tobramycin MIC (µg/mL)", fontsize=11, labelpad=10
)
# Add a light grid.

ax.grid(
    alpha=0.25,
    zorder=0,
)
fig.subplots_adjust(bottom=0.15)

# Add the colour scale for the number of isolates.

colour_bar = fig.colorbar(
    scatter,
    ax=ax,
    pad=0.02,
)

colour_bar.set_label(
    "Number of isolates"
)


# Create the statistical summary without a p-value.

statistics_text = (
    f"Spearman ρ\n"
    f"{results['spearman_rho']:.3f}\n\n"
    f"Bootstrap SE\n"
    f"{results['bootstrap_standard_error']:.3f}\n\n"
    f"95% BCa interval\n"
    f"[{results['bca_95_lower']:.3f}, "
    f"{results['bca_95_upper']:.3f}]\n\n"
    f"Sample size\n"
    f"n = {int(results['n'])} isolates"
)


# Place the statistical summary in the separate right-side panel.

information_panel.text(
    0.05,
    0.70,
    statistics_text,
    horizontalalignment="left",
    verticalalignment="center",
    fontsize=11,
    linespacing=1.2,
    bbox={
        "boxstyle": "round",
        "facecolor": "white",
        "edgecolor": "red",
        "linewidth": 1.5,
        "alpha": 0.95,
    },
)


# Explain the bubble labels on the right side of the figure.
fig.text(
    0.70,
    0.30,
    "Bubble size and labels show\n"
    "how many isolates share the\n"
    "same MIC–BPC combination.",
    horizontalalignment="left",
    verticalalignment="center",
    fontsize=9,
    color="dimgray"
)


# Save and close Figure 1.

fig.savefig(
    "figure1_MIC_BPC_association.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


# 5. Create Figure 2: bootstrap uncertainty.

fig, ax = plt.subplots(
    figsize=(10, 6)
)


# Plot the distribution of the bootstrap correlation values.

ax.hist(
    bootstrap_values,
    bins=35,
    color="lightblue",
    edgecolor="black",
    alpha=0.85,
)


# Add the observed Spearman correlation.

ax.axvline(
    results["spearman_rho"],
    color="red",
    linewidth=2,
    label=(
        f"Observed ρ = "
        f"{results['spearman_rho']:.3f}"
    ),
)


# Add the lower limit of the BCa uncertainty interval.

ax.axvline(
    results["bca_95_lower"],
    color="green",
    linestyle="--",
    linewidth=2,
    label="95% BCa interval",
)


# Add the upper limit of the BCa uncertainty interval.

ax.axvline(
    results["bca_95_upper"],
    color="green",
    linestyle="--",
    linewidth=2,
)


# Shade the region between the interval limits.

ax.axvspan(
    results["bca_95_lower"],
    results["bca_95_upper"],
    color="green",
    alpha=0.10,
)


# Add the Figure 2 title and axis labels.

ax.set_title(
    "Uncertainty in the Spearman Correlation",
    fontweight="bold",
    pad=15,
)

ax.set_xlabel(
    "Bootstrapped Spearman correlation (ρ)",
    labelpad=10,
)

ax.set_ylabel(
    "Frequency",
    labelpad=10,
)


# Add a light horizontal grid and the legend.

ax.grid(
    axis="y",
    alpha=0.25,
)

ax.legend(
    loc="upper left",
)


# Prevent labels from overlapping with the plot.

fig.tight_layout()


# Save and close Figure 2.

fig.savefig(
    "figure2_bootstrap_uncertainty.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


# 6. Confirm that both figures were created.

print("Both Matplotlib figures were created successfully.")

print(
    "Figure 1: figure1_MIC_BPC_association.png"
)

print(
    "Figure 2: figure2_bootstrap_uncertainty.png"
)
