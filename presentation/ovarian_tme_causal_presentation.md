---
marp: true
theme: default
paginate: true
title: Ovarian TME Causal Structure Learning
---

# Ovarian TME Causal Structure Learning

From single-cell ovarian cancer data to hypothesis-generating causal graphs

Base reference:
Zheng et al. Nature Cancer 2023

Project goal:
- Recast ovarian TME correlations as a constrained causal discovery problem
- Compare primary tumor, metastatic tumor, lymph node, ascites, and PBMC ecosystems
- Identify stable edges linking tumor burden, stromal context, myeloid states, and lymphocyte states

---

# Study Motivation

Why this analysis:
- The Zhang/Zheng ascites study described strong correlative remodeling across tissue compartments
- We want to test whether those relationships organize into a plausible directed structure
- The current analysis is best viewed as hypothesis generation, not definitive causal proof

Working biological questions:
- Does tumor-heavy tissue co-occur with suppressive immune states?
- Do dendritic-cell-rich myeloid states align with adaptive immune infiltration?
- Are cytotoxic programs separable from broader immune-infiltrated programs?

---

# Dataset Overview

Summary from the processed metadata:
- 14 patients
- 39 samples
- 5 tissue compartments
- Broad cell annotations used to engineer sample-level features

Tissue distribution:
- Primary Tumor: 61,302 cells
- Ascites: 53,499 cells
- Lymph Node: 40,974 cells
- PBMC: 35,501 cells
- Metastatic Tumor: 32,087 cells

Data source tables:
- [tab01_dataset_overview.csv](/Users/aingela/ovarian-tme-causal/tables/tab01_dataset_overview.csv)
- [tab03_cells_per_tissue.csv](/Users/aingela/ovarian-tme-causal/tables/tab03_cells_per_tissue.csv)

---

# Broad TME Composition

Most abundant broad cell classes:
- CD4+ T: 47,440 cells
- Epithelial cells: 44,805 cells
- Macrophage: 34,523 cells
- CD8+ T: 29,573 cells
- B: 18,057 cells
- Monocyte: 11,075 cells
- NK: 10,125 cells

Interpretation:
- The dataset contains both tumor-rich and immune-rich samples
- Macrophages, T cells, and epithelial cells are all abundant enough to drive higher-order TME modules

![Broad Composition](/Users/aingela/ovarian-tme-causal/figures/fig08_composition_per_tissue.png)

---

# UMAP Overview

This slide anchors the feature engineering to the single-cell structure.

What to emphasize:
- Broad cell states are visually well separated
- Tissue-specific composition shifts are visible across the atlas
- The sample-level causal features are derived from this cell-state structure

![UMAP broad cell types](/Users/aingela/ovarian-tme-causal/figures/fig03_umap_maintypes2.png)

---

# T Cell and Myeloid State Design

T-cell states used in the feature matrix:
- Cytotoxic
- Exhausted
- Regulatory
- Naive
- Effector-memory

Myeloid states used in the feature matrix:
- TAM immunosuppressive
- Dendritic cell
- Monocyte
- Macrophage other

These states were defined to map onto interpretable ovarian TME programs described in the literature.

Supporting plots:
- ![T cell features](/Users/aingela/ovarian-tme-causal/figures/fig13_tcell_cyto_exh_features.png)
- ![T cell marker heatmap](/Users/aingela/ovarian-tme-causal/figures/fig12_tcell_marker_heatmap.png)

---

# Feature Matrix for DAG Learning

Selected sample-level variables:
- `prop_CD8_T`
- `prop_CD4_T`
- `prop_Epithelial_cells`
- `prop_Fibroblast`
- `prop_NK`
- `prop_B`
- `prop_Endothelial_cells`
- `CD8_CD4_ratio`
- `Immune_Tumor_ratio`
- `Tcell_Myeloid_ratio`
- `Tstate_simple_exhausted`
- `Tstate_simple_cytotoxic`
- `Tstate_simple_regulatory`
- `Mstate_simple_TAM_immunosuppressive`
- `Mstate_simple_dendritic_cell`

Important constraint:
- 39 samples vs 15 variables gives an n/p ratio of about 2.6
- This means the graphs are exploratory and likely to contain redundancy-driven ambiguity

---

# Correlation Structure Reveals TME Modules

The selected features cluster into biologically coherent blocks:
- Structural/stromal: fibroblast, endothelial, CD8/CD4 ratio
- Tumor-suppressive: epithelial, TAM immunosuppressive, exhausted, regulatory
- APC/helper/infiltrated: B cells, dendritic cell, Tcell/Myeloid ratio, CD4, Immune/Tumor ratio
- Cytotoxic effector: CD8, NK, cytotoxic state

Take-home message:
- The data are not random; they organize into higher-order ecosystem programs
- This is encouraging for downstream DAG learning, but it also creates redundancy

![Feature correlations](/Users/aingela/ovarian-tme-causal/figures/fig_feature_correlations.png)

---

# Tissue-Stratified Biology

Examples of tissue-enriched patterns:
- Primary and metastatic tumors are enriched for epithelial burden and TAM-immunosuppressive states
- Lymph node samples are enriched for dendritic-cell and B-cell features
- Ascites and PBMC show stronger immune-skewed signals than tumor-heavy sites

Illustrative means by tissue:
- `Mstate_simple_TAM_immunosuppressive`: highest in primary and metastatic tumor
- `Mstate_simple_dendritic_cell`: highest in lymph node
- `Tstate_simple_exhausted`: highest in primary and metastatic tumor
- `prop_Epithelial_cells`: highest in primary tumor

Data source:
- [tab_features_by_tissue.csv](/Users/aingela/ovarian-tme-causal/tables/tab_features_by_tissue.csv)

---

# Causal Discovery Strategy

Methods run on the selected feature matrix:
- PC algorithm: constraint-based
- FCI: allows latent confounding
- Hill-climbing: score-based
- TABU: score-based with broader search
- Bootstrap averaging: stability assessment

Interpretation framework:
- Preliminary graphs are permissive and used for signal discovery
- Stringent graphs are used to look for reproducible structure
- Bootstrap support is the main filter for robustness

Method timing was light locally, so the bottleneck is inference quality, not compute:
- PC: 0.036 sec
- FCI: 0.019 sec
- HC: 0.006 sec
- TABU: 0.006 sec
- Bootstrap 50 iters: 1.83 sec

Generated causal figure files:
- [fig_01_pc_preliminary.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_01_pc_preliminary.pdf)
- [fig_05_pc_stringent.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_05_pc_stringent.pdf)
- [fig_07_bootstrap_strength_heatmap.png](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_07_bootstrap_strength_heatmap.png)
- [fig_08_bootstrap_stability_barplot.png](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_08_bootstrap_stability_barplot.png)
- [fig_09_final_causal_graph.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_09_final_causal_graph.pdf)
- [fig_10_exploratory_causal_graph.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_10_exploratory_causal_graph.pdf)

---

# Reading The Causal Figures

How to present the final causal output folder:
- `fig_01` and `fig_05` show PC-based graph structure before and after stringent filtering
- `fig_07` shows edge strength across all ordered pairs under bootstrap
- `fig_08` ranks the strongest bootstrap-supported edges
- `fig_09` is the final bootstrap-averaged causal graph at the 50% threshold
- `fig_10` is the more permissive exploratory graph at the 30% threshold

Recommended framing:
- Use `fig_10` to explain overall module structure
- Use `fig_09` to report the short list of stable, reproducible relationships
- Use `fig_07` and `fig_08` to explain why some arrows remain ambiguous

---

# Bootstrap Edge Strength Heatmap

File:
- [fig_07_bootstrap_strength_heatmap.png](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_07_bootstrap_strength_heatmap.png)

What this figure shows:
- A dense matrix of directed edge support across bootstrap runs
- Hot cells mark edges that repeatedly reappear
- Reciprocal hot cells indicate stable association but unstable direction

What this file says about the final DAG:
- The dominant signals are concentrated in a few feature pairs, not spread evenly across the graph
- Several of the strongest blocks are reciprocal, which is consistent with redundancy among selected variables
- The final graph should therefore be interpreted as a sparse summary of stronger modules, not as a complete mechanistic map

![Bootstrap edge heatmap](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_07_bootstrap_strength_heatmap.png)

---

# Stable Bootstrap Edges

Top bootstrap-supported relationships at the 50% threshold:
- `prop_CD4_T <-> Mstate_simple_TAM_immunosuppressive` with strength 0.999
- `Tcell_Myeloid_ratio <-> Mstate_simple_dendritic_cell` with strength 0.973
- `prop_Fibroblast <-> prop_Endothelial_cells` with strength 0.933
- `prop_NK <-> Tstate_simple_cytotoxic` with strength 0.927
- `prop_Endothelial_cells <-> CD8_CD4_ratio` with strength 0.846

How to read this:
- The relationships are stable
- The directions are often not stable
- Bidirectional support usually means the pair is biologically linked but not orientable with this sample size

![Bootstrap stability](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_08_bootstrap_stability_barplot.png)

---

# High-Confidence Consensus Edges

Edges supported by both stringent PC and bootstrap:
- `Mstate_simple_dendritic_cell -> Tcell_Myeloid_ratio`
- `Mstate_simple_dendritic_cell -> prop_B`
- `prop_CD4_T -> Mstate_simple_TAM_immunosuppressive`
- `prop_Endothelial_cells -> CD8_CD4_ratio`
- `prop_Epithelial_cells -> Tstate_simple_regulatory`
- `prop_Fibroblast -> prop_Endothelial_cells`
- `prop_NK -> Tstate_simple_cytotoxic`

Biologic interpretation:
- A dendritic/B/helper axis is separable from the tumor-suppressive axis
- Tumor-rich samples align with regulatory/suppressive T-cell biology
- Stromal and vascular features travel together
- NK abundance and cytotoxic state behave as a shared effector module

Data source:
- [final_edge_consensus.csv](/Users/aingela/ovarian-tme-causal/results/causal/final_edge_consensus.csv)

---

# Final DAG Result

Primary final graph file:
- [fig_09_final_causal_graph.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_09_final_causal_graph.pdf)

This final DAG is represented by the stable consensus edges:
- `prop_Fibroblast -> prop_Endothelial_cells`
- `prop_Endothelial_cells -> CD8_CD4_ratio`
- `prop_Epithelial_cells -> Tstate_simple_regulatory`
- `prop_CD4_T -> Mstate_simple_TAM_immunosuppressive`
- `Mstate_simple_dendritic_cell -> prop_B`
- `Mstate_simple_dendritic_cell -> Tcell_Myeloid_ratio`
- `prop_NK -> Tstate_simple_cytotoxic`

Interpretation of the final structure:
- A structural tissue axis links fibroblast and endothelial abundance
- A tumor-dominant axis links epithelial burden with regulatory T-cell enrichment
- A dendritic-cell-centered immune organization axis links dendritic signal with B cells and T cell to myeloid balance
- A cytotoxic effector axis links NK abundance with the cytotoxic T-cell state feature

How to present it:
- Treat this as the most conservative graph in the project
- Emphasize that these are the edges surviving both bootstrap support and stringent PC support
- Note that some biologically plausible pairs did not survive with stable direction, which is expected at this sample size

---

# Exploratory vs Final Graphs

Exploratory graph:
- Better for visualizing the broader ecosystem structure
- Includes weaker and more redundant edges

Final graph:
- Sparse and more defensible
- Best used for reporting the small set of reproducible hypotheses

Recommended interpretation:
- Report modules and stable edges
- Avoid claiming strong directionality for pairs with near-symmetric support
- Use the final graph as a shortlist for follow-up validation

Files:
- [fig_10_exploratory_causal_graph.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_10_exploratory_causal_graph.pdf)
- [fig_09_final_causal_graph.pdf](/Users/aingela/ovarian-tme-causal/results/causal/figures/fig_09_final_causal_graph.pdf)

Exploratory graph highlights from preliminary consensus:
- `Mstate_simple_TAM_immunosuppressive -> Tstate_simple_exhausted`
- `Mstate_simple_dendritic_cell -> prop_B`
- `prop_Endothelial_cells -> CD8_CD4_ratio`
- `prop_Fibroblast -> prop_Endothelial_cells`
- `prop_NK -> prop_CD8_T`

Final graph highlights from stringent consensus:
- Keeps the dendritic/B-cell axis
- Keeps the fibroblast/endothelial structural axis
- Keeps the epithelial/regulatory and NK/cytotoxic links
- Drops several weaker suppressive-module arrows or leaves them directionally unresolved

---

# Main Biological Takeaways

What the current results support:
- Tumor-heavy compartments co-localize with suppressive immune programs
- Dendritic-cell-rich myeloid programs align with B-cell and immune-infiltrated structure
- Cytotoxic activity forms a distinct module centered on NK and CD8-like effector biology
- Stromal and endothelial features behave as a structural tissue program

What remains uncertain:
- Exact arrow direction for several highly correlated pairs
- Whether some edges reflect composition rather than mechanism
- Whether patient-level or site-level latent confounders are still driving parts of the graph

---

# Limitations

Key caveats:
- Small sample size relative to variable count
- Proportion data are compositional
- Several selected features are strongly correlated by construction
- The DAG is cross-sectional and cannot establish temporal causality

Best framing:
- This is a causal discovery screen that prioritizes mechanistically plausible relationships
- The strongest value is in refining which correlative observations from Zheng et al. deserve focused validation

---

# Recommended Next Steps

Analysis improvements:
- Reduce redundancy by using one representative per module
- Add domain-knowledge constraints to the DAG search
- Test module-level DAGs rather than only raw features
- Include tissue site as an explicit upstream context node

Biological follow-up:
- Validate tumor-regulatory/TAM links
- Validate dendritic-B-helper module across tissue compartments
- Test whether cytotoxic and APC/helper modules are separable in external datasets

Project files:
- [results/causal](/Users/aingela/ovarian-tme-causal/results/causal)
- [figures](/Users/aingela/ovarian-tme-causal/figures)
- [tables](/Users/aingela/ovarian-tme-causal/tables)
