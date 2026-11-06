# Checkpoint 1
To run this:

- Open `project_checkpoint_notebook.ipynb`, select the installed Python environment, and run all cells.
- `ONLY_22 = True` runs all three questions on chromosome 22 to save time. Set it to `False` to reproduce the paper’s chromosome 1–22 eQTL and eGene counts. The EUR–YRI comparison always uses chromosome 22. Data for both modes is included.

## If you want to run from a terminal

Open a terminal in this folder.

macOS ARM:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
chmod +x tools/plink2
python run_checkpoint.py
```

Windows 64-bit (PowerShell):

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_checkpoint.py
```

The correct bundled PLINK executable is selected automatically.

## Data and results

- All data needed to run the notebook is included in `checkpoint_data/`. No additional datasets need to be downloaded.
- Source details are in the included `sources.json` files.
- Results appear below the notebook cells and in `results/`.
- The terminal runner saves `results/executed_checkpoint.ipynb`.
