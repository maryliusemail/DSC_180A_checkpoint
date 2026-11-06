# Checkpoint 1
To run this:

- Open `project_checkpoint_notebook.ipynb`, select the installed Python environment, and run all cells.
- `TA_GRADING = True` runs on chromosome 22 by default. If it's set to `False`, it will run on all chromosomes (1–22), which the data files are too large to be uploaded to GitHub.

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
