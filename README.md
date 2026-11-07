# Checkpoint 1

<<<<<<< HEAD
## Setup
=======
Open `project_checkpoint_notebook.ipynb` and run all cells. All of the analysis codes are in the notebook.
The tags are:
- `ONLY_22 = True`: means chromosome 22 only.
- `ONLY_22 = False`: means chromosomes 1–22. The EUR–YRI comparison always uses chromosome 22.
>>>>>>> 7ab10dc985e2a483d4b074e5220bf846d36764ba

Use Python 3.12. From this repository's folder, install the packages:

```sh
python -m pip install -r requirements.txt
```

PLINK 2 is included in `tools/`: `plink2` for Apple Silicon Macs and `plink2.exe` for Windows. 

On Mac/Linux, run:

```sh
chmod +x tools/plink2
```

## Data

All the required data inputs are in `checkpoint_data/`. If you want to see the sources, they are in `checkpoint_data/sources.json`.

## Run

From the repository's folder:

```sh
python -m jupyter lab project_checkpoint_notebook.ipynb
```

Open the notebook and choose Run All Cells.

- ONLY_22 = False: means that it will run the analysis on chromosomes 1–22, which is what is described in the report.
- ONLY_22 = True: means that it will run the analysis on only chromosome 22. This will save time.
- The EUR–YRI comparison always uses chromosome 22. 

The results appear in the notebook as tables or words, and also can be seen in  `results/`. If you rerun, it will replace the past results.
