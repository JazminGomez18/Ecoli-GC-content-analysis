# GC Content Analysis in *E. coli*

This project analyzes the GC content (%GC) of *E. coli* genes and compares it with random expectations based on the genome composition.  
The goal is to determine how likely each gene’s GC% is to occur by chance, using an empirical randomization test.

---

## Method summary

- The *E. coli* genome sequence (`e_coli.fasta`) and the gene sequences (`genes_ecoli.txt`) are loaded.  
- For each gene:
  - The observed GC% is calculated.
  - 500 random fragments of the same length are extracted from the genome.
  - The mean GC% of these fragments is used as the expected value.
  - A two-sided empirical p-value is calculated by comparing the observed and expected GC%.
- Results are saved to a table and visualized in Seaborn plots.

---

## Outputs

| File | Description |
|------|--------------|
| `GC_per_gene.tsv` | Table with `Gene`, `Length`, `GC_observed`, `GC_random_mean`, and `p_value`. |
| `GC_distribution_genes.png` | Histogram showing the distribution of GC% across all genes. |
| `GC_vs_random_simple.png` | Scatterplot comparing observed vs. expected GC% for each gene. |

---
