# Credits, starting materials and third-party notices

**This effort is supported by Qollab & IonQ**

## Research methods and prior project

- **Wataru Onodera, Nobuyuki Hara, Shiho Aoki, Toru Asahi and Naoya Sawamura**.
  *Phylogenetic tree reconstruction via graph cut presented using a
  quantum-inspired computer*. Molecular Phylogenetics and Evolution 178 (2023),
  107636; online October 2022.
  [Paper](https://doi.org/10.1016/j.ympev.2022.107636).
  Credit for NMcutDA, recursive normalized-cut methodology, penalty
  reinforcement and the model-NJ comparison motivation. Their solver was
  Fujitsu's **Digital Annealer, a quantum-inspired classical machine**.
- **Nicolás Fernández-Otero, Tomás F. Pena and Juan C. Pichel**.
  *Reconstruction of phylogenetic trees via graph-splitting using quantum
  computing*. The Journal of Supercomputing 82 (2026), 324.
  [Paper](https://doi.org/10.1007/s11227-026-08465-x).
  Credit for the earlier QA/QAOA phylogeny approach that motivated this project.
  [Nicofero/PhyloTree_ReConst](https://github.com/Nicofero/PhyloTree_ReConst)
  supplied the inspected implementation and selected alignment/reference files.
  Local checkout inspected: `2fba9f1904a13fa41246af9cef25574f367f0239`.
- **Edward Farhi, Jeffrey Goldstone and Sam Gutmann**. *A Quantum Approximate
  Optimization Algorithm* (2014). [arXiv:1411.4028](https://arxiv.org/abs/1411.4028).
  Credit for QAOA.
- **Qollab tutorial authors**: the injected `backend`, top-level circuit execution,
  job polling and tutorial examples informed platform integration. The final
  runner does not import the tutorial optimizer's reported probabilities.

Bibliographic details above were checked against the locally supplied papers
and upstream repository. Research PDFs are referenced rather than redistributed.
No paper text or figure is presented as an original project contribution.

## What this hackathon added

The starting materials were the user's prior planning/handoff, initial
`qphylo_core.py` proposal, the forked repository, research papers and Qollab
tutorial examples. QPhylo does not claim to invent QAOA or NMcut phylogeny.

The hackathon work implements a reusable Qiskit sampled-cut adapter; explicit
bit-order/energy checks; local induced-graph normalization; bounded reinforcement;
coordinate/COBYLA pilot comparison; an audited BLASTP affinity bridge; model-NJ
and bootstrap comparison; and a self-contained runner that respects Qollab's
backend job limits. It records counts, failures, decisions and hashes so that
optimization error can be distinguished from an affinity-method discrepancy.
The new modules are a separate implementation informed by the cited methods and
inspected upstream work. The nested fork remains separate and was not rewritten.

Noah directed the research question, reviewed methods, selected settings and ran
the Qollab experiments. OpenAI Codex assisted with code, testing, analysis and
documentation under that direction. This AI-assisted development is disclosed;
saved Qollab measurements are experimental evidence, while mocked tests and exact
calculations are explicitly separate.

The optional `project_banner_v2.png` header was created with the built-in OpenAI
image generation tool. It is conceptual protein/quantum artwork, not a measured
tree, validated molecular structure or actual circuit diagram. No stock image or
paper figure was supplied. The scientific results figures are separate original
plots of saved project data. The generation prompt is retained in the source
release's `docs/HEADER_ART.md`.

## Selected data and reference

The copied input is `PhyloTree_ReConst/alignments/8/assembled_sequences_1822.fasta`
(SHA-256 `9c80c2e7af4fe29dc6edda1d291204b27f22d81168169189c40f8f7d805ff24b`).
The supplied comparison tree is `PhyloTree_ReConst/trees/8/tree_best_1822.newick`
(SHA-256 `73d9dfe81c42a849ad5e6cfcbcdc4fc44d7695dccd86478e63636560f08bbd3c`).
Both are retained unchanged with the upstream repository's MIT notice below.

The FO paper identifies [RAxMLGrove](https://github.com/angtft/RAxMLGrove) as its
sequence/tree source. This does not establish the generating history, original
alignment source, simulation parameters or source-specific terms of this exact
1822 pair. File names and notebook labels alone are insufficient. We call this
tree the **provided reference**, not confirmed evolutionary ground truth. The
repository MIT notice is preserved; it is not asserted to settle every earlier
dataset's rights. Derived affinities and NJ results record the selected input
hash and preprocessing. Original article/data ownership remains with its authors.

RAxMLGrove credits **D. Höhler, W. Pfeiffer, V. Ioannidis, H. Stockinger and
A. Stamatakis (2022)**, *RAxML Grove: an empirical phylogenetic tree database*,
Bioinformatics 38(6):1741–1742. [Paper](https://doi.org/10.1093/bioinformatics/btab863).
Its [repository notice](https://github.com/angtft/RAxMLGrove#readme) specifies the
[Open Database License 1.0](https://opendatacommons.org/licenses/odbl/1-0/) for the
database and [Database Contents License 1.0](https://opendatacommons.org/licenses/dbcl/1-0/)
for individual contents. These source notices are retained here; upstream data
are not relicensed as QPhylo's original MIT software. We have not independently
mapped the selected files back to a particular Grove entry or inferred their
generating-tree status.

## Software dependencies

Dependencies are used, not vendored or relicensed by QPhylo. Their own licenses
apply: [Qiskit](https://github.com/Qiskit/qiskit) (Apache-2.0),
[SciPy](https://github.com/scipy/scipy) (BSD-3-Clause),
[NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/),
[R](https://www.r-project.org/) (GPL),
[ape](https://cran.r-project.org/package=ape) (GPL-2/GPL-3), and
[phangorn](https://cran.r-project.org/package=phangorn) (GPL-2 or later).
Qiskit/provider execution is supplied by Qollab; IonQ simulator measurements are
not hardware runs. Python, Conda, Jupyter and their maintainers are also credited.

## Preserved PhyloTree_ReConst notice

The following notice applies to included upstream material and is reproduced
verbatim from its `LICENSE`. QPhylo's MIT license does not replace it.

```text
MIT License

Copyright (c) 2024 nicofero01@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
