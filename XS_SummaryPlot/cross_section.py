import matplotlib.pyplot as plt
import numpy as np

# Central values in fb (add the new σ=3.11)
y = np.array([0.56, 1.11, 2.92, 2.73, 2.89, 3.11])

# Statistical uncertainties
stat_up = np.array([0.67, 0.41, 0.48, 0.22, 0.53, 0.22])
stat_down = np.array([0.44, 0.35, 0.44, 0.22, 0.49, 0.22])

# Systematic uncertainties
syst_up = np.array([0.21, 0.14, 0.28, 0.15, 0.29, 0.14])
syst_down = np.array([0.06, 0.10, 0.24, 0.15, 0.21, 0.12])

# Compute total uncertainties in quadrature
total_up = np.sqrt(stat_up**2 + syst_up**2)
total_down = np.sqrt(stat_down**2 + syst_down**2)

# Compute relative uncertainties (%) 
total_rel_up = total_up / y * 100
total_rel_down = total_down / y * 100
syst_rel_up = syst_up / y * 100
syst_rel_down = syst_down / y * 100

# X positions
x = np.arange(6)

# Labels with luminosity or measurement info
labels = [
    "7 TeV\n(5.1 fb$^{-1}$)",
    "8 TeV\n(19.7 fb$^{-1}$)",
    "13 TeV\n(35.9 fb$^{-1}$)",
    "13 TeV\n(Run2 UL)",
    "13.6 TeV\n(34.7 fb$^{-1}$)",
    "13.6 TeV\n(171 fb$^{-1}$)"
]

# Plot
fig, ax = plt.subplots(figsize=(12,5))

# Horizontal reference line
ax.axhline(0, linestyle='--', color='gray', linewidth=1)

# Black error bars: total relative uncertainty
ax.errorbar(
    x, [0]*6,
    yerr=[total_rel_down, total_rel_up],
    fmt='o',
    color='black',
    capsize=5,
    elinewidth=2,
    markersize=6,
    label="Total uncertainty"
)

# Red error bars: systematic only
ax.errorbar(
    x, [0]*6,
    yerr=[syst_rel_down, syst_rel_up],
    fmt='none',
    ecolor='red',
    elinewidth=6,
    capsize=0,
    alpha=0.9,
    label="Systematic only"
)

# Make all spines black and visible (full border)
for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')
    spine.set_visible(True)

# Axes formatting
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=14)
ax.set_ylabel("Relative uncertainty (%)", fontsize=16)

# Center y-axis around 0 with symmetric range
ymax = max(max(total_rel_up), max(total_rel_down))
ax.set_ylim(-ymax*1.1, ymax*1.1)  # 10% margin

# Add legend
ax.legend(fontsize=12, loc='upper right')

plt.tight_layout()
plt.savefig("uncertainty_plot_cms_updated.png", dpi=300, bbox_inches='tight')
plt.show()