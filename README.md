# Project checkpoint — TA check

Open `project check point.ipynb`, leave **`TA_GRADING = True`**, and run all cells.

The TA check uses only the inputs in this folder and tests **chromosome 22**:

- **344** matched individuals.
- **280,105** SNP–gene tests.
- **1,317** significant pairs after chromosome-wide Bonferroni correction.
- **72 eGenes** out of **633** tested genes.

## Contents

```text
checkpoint/
  project check point.ipynb
  requirements.txt
  checkpoint_data/
    LDREF/          # Copied chromosome 22 BED/BIM/FAM files
    expression/     # All chromosome 22 expression rows, original sample columns
    annotations/    # Copied Week 4 GRCh37 gene records
    sources.json    # Source paths, subset description, and SHA-256 checksums
  tools/plink2      # Copied executable, macOS ARM build
  results/ta_chr22/ # Generated results
```

Original project data remain under `../data/`. The expression subset preserves the original text values for chromosome 22, rather than re-normalizing them. No personal genotype data are included.

## Run

With Python 3.12 and this folder as the working directory:

```sh
python -m pip install -r requirements.txt
python -m jupyter lab "project check point.ipynb"
```

Select the Python environment containing those dependencies, then restart the kernel and run all cells. The default TA check does not require the original project data folder or GENCODE download. Results are written under `results/ta_chr22/`; rerunning replaces them.

The supplied PLINK executable targets macOS ARM. For another platform, place a compatible PLINK 2 executable at `tools/plink2` and ensure it is executable.

## Full analysis

Set `TA_GRADING = False` and run all cells to test chromosomes 1–22. This mode requires the original project's `../data/` folder, including all 22 LDREF panels, the full expression matrix, and GENCODE v12. Its outputs remain in `../results/project_checkpoint/`. The small TA data folder alone supports chromosome 22 only.

The model reproduces Week 4 unadjusted regression; it does not apply PCA correction. Chromosome-wide and scan-wide Bonferroni thresholds coincide when only chromosome 22 is selected.

Verified on November 2, 2026: all five cells passed from this folder and from an isolated copy without the original project data. Both runs reproduced the expected chromosome 22 results.

A separate local notebook remains at `../notebooks/project check point.ipynb`. It uses the original project data and tools, so moving this submission folder does not break the local copy.
