import os
import subprocess
import tempfile

# Snakemake provides a global `snakemake` object when using the script: directive
# Expected rule interface:
# - input: results (list of paths)
# - params: r_script (path to the R script to run)
# - output: single PDF path
# - conda: r_plotting.yaml (provides Rscript and Python)

def _touch(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as fh:
        pass

def main():
    # snakemake.input may be a sequence; when named, attribute access works
    results_attr = getattr(snakemake.input, "results", None)
    results = list(results_attr) if results_attr is not None else list(snakemake.input)
    out_pdf = snakemake.output[0]
    r_script = getattr(snakemake.params, "r_script", None)

    if not results or len(results) <= 1:
        _touch(out_pdf)
        print("Not enough input files, writing placeholder and skipping.")
        return

    if not r_script or not os.path.exists(r_script):
        raise FileNotFoundError(f"R script not found: {r_script}")

    # Write results list to a temporary file to avoid ARG_MAX issues
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".list") as handle:
        handle.write("\n".join(map(str, results)))
        list_path = handle.name

    cmd = [
        "Rscript",
        r_script,
        "--results-file",
        list_path,
        "--output",
        out_pdf,
    ]

    try:
        subprocess.check_call(cmd)
    finally:
        try:
            os.remove(list_path)
        except OSError:
            pass

if __name__ == "__main__":
    main()
