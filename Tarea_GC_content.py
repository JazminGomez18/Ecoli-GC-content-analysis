# Import libraries
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random

# Load genome sequence
genome = ""
with open('../secuencias/e_coli.fasta') as file:
    next(file)
    for line in file:
        genome += line.strip()

# Load gene sequences
gene_dict = {}
with open('../secuencias/genes_ecoli.txt') as file:
    for line in file:
        if line.startswith(">"):
            temp = line.split()[1][6:-1]
            gene = temp
            gene_dict[gene] = ""
        else:
            gene_dict[gene] += line.strip()

# Define functions
def GC_content(sequence):
    count_A = sequence.count('A')
    count_G = sequence.count('G')
    count_C = sequence.count('C')
    count_T = sequence.count('T')
    GC_content = (count_G + count_C) / len(sequence)
    GC_content = round(100 * GC_content, 2)
    return GC_content
    pass


def gc_random_stats(sequence, genome_sequence, n_iter=1000):

    L = len(sequence)
    genome_len = len(genome_sequence)
    gc_list = []

    for _ in range(n_iter):
        start = random.randint(0, genome_len - L)
        fragment = genome_sequence[start:start + L]
        gc_value = GC_content(fragment)
        gc_list.append(gc_value)

    gc_series = pd.Series(gc_list)
    gc_mean = round(gc_series.mean(), 2)
    gc_obs = GC_content(sequence)

    # two-sided empirical p-value
    diff_obs = abs(gc_obs - gc_mean)
    diff_null = [abs(x - gc_mean) for x in gc_list]
    p_value = sum(d >= diff_obs for d in diff_null) / len(diff_null)

    return gc_obs, gc_mean, round(p_value, 5)

# 4. Compute GC stats per gene
results = []

for gene, seq in gene_dict.items():
    gc_obs, gc_mean, p_value = gc_random_stats(seq, genome, n_iter=500)
    results.append({
        "Gene": gene,
        "Length": len(seq),
        "GC_observed": gc_obs,
        "GC_random_mean": gc_mean,
        "p_value": p_value
    })

# 5. Save results and plot
df = pd.DataFrame(results)
df.to_csv("GC_per_gene.tsv", sep="\t", index=False)

plt.figure(figsize=(8,5))
sns.histplot(data=df, x="GC_observed", bins=40, kde=True, color="darkgreen", alpha=0.6)
plt.axvline(df["GC_observed"].mean(), color="black", linestyle="--", linewidth=2, label="Mean GC")
plt.title("Distribution of GC% across E. coli genes")
plt.xlabel("%GC")
plt.ylabel("Frequency")
plt.legend()
plt.tight_layout()
plt.savefig("GC_distribution_genes.svg", dpi=600)
plt.show()

plt.figure(figsize=(7,6))

sns.scatterplot(
    data=df,
    x="GC_random_mean",
    y="GC_observed",
    s=40,
    color="darkgreen",
    alpha=0.7
)

# Línea diagonal de referencia (observado = esperado)
plt.plot([30,70], [30,70], color="gray", linestyle="--", linewidth=1)

plt.title("Observed vs Random GC% per gene")
plt.xlabel("Expected GC% (random fragments)")
plt.ylabel("Observed GC%")

plt.tight_layout()
plt.savefig("GC_vs_random_simple.svg", dpi=600)
plt.xlim(40, 60)
plt.show()