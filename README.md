# Project checkpoint

All three questions use **chromosome 22**:

- How many significant SNP–gene pairs?
- How many eGenes (genes with at least one significant pair)?
- How many EUR–YRI differences in SNP effects?

## Run

All inputs are included in `checkpoint_data/`. Keep the folder structure intact.
Use Python 3.12. On macOS ARM:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
chmod +x tools/plink2
python run_checkpoint.py
```

On Windows 64-bit (PowerShell):

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_checkpoint.py
```

Run these commands from this folder. The notebook selects `tools/plink2.exe` on
Windows and `tools/plink2` on macOS. Both executables are included. Linux requires
a compatible [PLINK 2 build](https://www.cog-genomics.org/plink/2.0/).
The full analysis has been tested on macOS ARM; Windows execution has not been tested.

To run interactively:

```sh
python -m jupyter lab "project_checkpoint_notebook.ipynb"
```

Select the environment’s kernel and run all cells. Results are saved in
`results/ta_chr22_pca/` and `results/population_chr22_pca/`.
The command-line runner also saves `results/executed_checkpoint.ipynb`.
Rerunning overwrites generated results. A run takes about 1–3 minutes on the tested machine.

## Data and methods

- `checkpoint_data/LDREF/`: European chromosome 22 genotypes; 344 samples match expression.
- `checkpoint_data/expression/`: GEUVADIS expression for 633 chromosome 22 genes.
- `checkpoint_data/annotations/`: Week 4 GRCh37 gene boundaries; test windows extend 500 kb on either side.
- `checkpoint_data/covariates/eur_loco/`: five European genotype PCs computed without chromosome 22.
- `checkpoint_data/populations/chr22_EUR_YRI/`: individual genotypes and labels for 358 EUR and 87 YRI samples, with 65,266 shared common SNPs.

The course supplied LDREF and expression data. Population genotypes come from
[1000 Genomes Phase 3](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
Source details, filtering commands, and checksums are recorded in the `sources.json`
files under `checkpoint_data/`. The full VCF is not needed to run this folder.
`scripts/prepare_eur_loco_pcs.py` reproduces European PCs from the full source project;
the included scores are sufficient for this checkpoint.

Association tests adjust for five PCs. EUR–YRI comparisons test the SNP–population
interaction with HC3 standard errors. Their PCs use chromosome 22 markers outside
each gene’s test window, so adjustment is limited. All counts use Bonferroni correction
at 0.05. Correlated SNPs can produce multiple significant pairs for one signal.

## Results

| Analysis | Tests | Significant |
|---|---:|---:|
| SNP–gene pairs | 280,105 | 1,302 |
| eGenes | 633 genes | 74 |
| EUR–YRI interactions | 1,137,206 | 341 pairs across 21 genes |

Earlier results without PCA remain in `results/ta_chr22/` and
`results/population_chr22/`. Changes in discovery counts alone do not establish
that PCA reduced false positives.
