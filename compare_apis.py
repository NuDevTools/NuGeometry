import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# --------------------------
# Load data
# --------------------------

def load_file(path):
    df = pd.read_csv(
        path,
        header=None,
        sep=' ',
        names=["x", "y", "z", "material", "energy"]
    )
    return df

api1 = load_file("prob_test.txt")
api2 = load_file("prob_test2.txt")

# --------------------------
# Histogram comparison
# --------------------------

def plot_histograms(df1, df2, column):
    plt.figure(figsize=(8,5))
    sns.histplot(df1[column], stat="density", bins=50, alpha=0.5, label="API 1")
    sns.histplot(df2[column], stat="density", bins=50, alpha=0.5, label="API 2")
    plt.title(f"Comparison of {column}")
    plt.legend()
    plt.tight_layout()
    plt.show()

numeric_cols = ["x", "y", "z", "energy"]

for col in numeric_cols:
    plot_histograms(api1, api2, col)

# --------------------------
# Statistical comparisons
# --------------------------

def compare_numeric(df1, df2, column):
    print(f"\n===== {column} =====")
    
    data1 = df1[column]
    data2 = df2[column]
    
    # Means and std
    print("Mean ± Std:")
    print(f"API 1: {data1.mean():.5f} ± {data1.std():.5f}")
    print(f"API 2: {data2.mean():.5f} ± {data2.std():.5f}")
    
    # KS test (distribution comparison)
    ks_stat, ks_p = stats.ks_2samp(data1, data2)
    print(f"KS test: statistic={ks_stat:.5f}, p-value={ks_p:.5e}")
    
    # Welch t-test (difference in means)
    t_stat, t_p = stats.ttest_ind(data1, data2, equal_var=False)
    print(f"Welch t-test: statistic={t_stat:.5f}, p-value={t_p:.5e}")
    
    # Levene test (variance comparison)
    lev_stat, lev_p = stats.levene(data1, data2)
    print(f"Levene test: statistic={lev_stat:.5f}, p-value={lev_p:.5e}")

for col in numeric_cols:
    compare_numeric(api1, api2, col)

# --------------------------
# Material comparison (categorical)
# --------------------------

print("\n===== Material Composition =====")

mat1 = api1["material"].value_counts()
mat2 = api2["material"].value_counts()

materials = sorted(set(mat1.index).union(set(mat2.index)))

counts1 = np.array([mat1.get(m, 0) for m in materials])
counts2 = np.array([mat2.get(m, 0) for m in materials])

contingency = np.vstack([counts1, counts2])

chi2, p, dof, expected = stats.chi2_contingency(contingency)

print("Material counts:")
for i, m in enumerate(materials):
    print(f"{m}: api1={counts1[i]}, api2={counts2[i]}")

print(f"\nChi-square test: chi2={chi2:.5f}, p-value={p:.5e}")

