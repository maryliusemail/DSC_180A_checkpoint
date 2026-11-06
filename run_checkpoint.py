"""Execute every code cell in the TA notebook without opening Jupyter."""
from pathlib import Path
import argparse
import base64
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/executed_checkpoint.ipynb",
                        help="Executed notebook path, relative to this folder")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    os.chdir(root)
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("IPYTHONDIR", str(root / "results/.ipython"))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import nbformat
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output

    notebook = nbformat.read(root / "project_checkpoint_notebook.ipynb", as_version=4)
    shell = InteractiveShell.instance()
    print("Runing the notebook", flush=True)
    cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    for count, cell in enumerate(cells, 1):
        cell.outputs = []
        cell.execution_count = count
        with capture_output() as captured:
            result = shell.run_cell(cell.source, store_history=True)
        for name, value in (("stdout", captured.stdout), ("stderr", captured.stderr)):
            if value:
                cell.outputs.append(nbformat.v4.new_output("stream", name=name, text=value))
        for rich in captured.outputs:
            data = {key: base64.b64encode(value).decode() if isinstance(value, bytes) else value
                    for key, value in rich.data.items()}
            cell.outputs.append(nbformat.v4.new_output("display_data", data=data, metadata=rich.metadata))
        nbformat.write(notebook, output)
        error = result.error_before_exec or result.error_in_exec
        if error:
            print(captured.stdout + captured.stderr, end="", flush=True)
            raise RuntimeError(f"Code cell {count} failed: {error}") from error
        if count == 5:
            print("eQTL and eGene counts complete.", flush=True)
        elif count == 8:
            print("EUR–YRI comparison complete.", flush=True)
    print(f"Done. Check the notebook: {output}", flush=True)


if __name__ == "__main__":
    main()
