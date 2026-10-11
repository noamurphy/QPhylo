# Protein1822: QAOA, exact NMcut and model-based NJ

NJ used WAG+Gamma (5 categories, shape0.5), 100 site bootstraps, seed42.
Input: the same 8 proteins / 389 aligned columns used for the QAOA graph.

| Method 1 | Method 2 | Unrooted RF | Normalized RF |
| --- | --- | --- | --- |
| QAOA | Exact_NMcut | 0 | 0 |
| QAOA | Model_NJ | 2 | 0.2 |
| QAOA | Provided_reference | 2 | 0.2 |
| Exact_NMcut | Model_NJ | 2 | 0.2 |
| Exact_NMcut | Provided_reference | 2 | 0.2 |
| Model_NJ | Provided_reference | 0 | 0 |

QAOA reproduces the exact affinity-tree topology on this run. Its disagreement with NJ is also present in exact NMcut, so it is not a QAOA topology approximation error on this input.

RF0 means identical unrooted splits. Reference agreement is conditional on the supplied reference; its generating-tree provenance is unverified.

## NJ site-bootstrap support

| Split sides | Support (%) |
| --- | --- |
| taxon2,taxon3 / taxon1,taxon4,taxon5,taxon6,taxon7,taxon8 | 100 |
| taxon4,taxon5 / taxon1,taxon2,taxon3,taxon6,taxon7,taxon8 | 100 |
| taxon1,taxon6,taxon7,taxon8 / taxon2,taxon3,taxon4,taxon5 | 89 |
| taxon1,taxon8 / taxon2,taxon3,taxon4,taxon5,taxon6,taxon7 | 94 |
| taxon6,taxon7 / taxon1,taxon2,taxon3,taxon4,taxon5,taxon8 | 100 |

NJ produced 2 negative branch lengths. These are retained in the Newick files; topology comparisons and the figure ignore branch lengths.

NJ bootstrap percentages describe site-resampling support, not probabilities that clades are correct.

The figure compares topology only: branch lengths are hidden and the last taxon is a display anchor, not a biological outgroup.

This reuses the provided alignment and records current R/library versions; it is not a complete reproduction of Onodera's alignment/tool versions.

![Tree comparison](tree_comparison.png)
