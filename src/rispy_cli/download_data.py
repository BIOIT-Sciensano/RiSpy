import os
import re
import shutil
import stat
import subprocess
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

DATA_URL = ("https://github.com/BIOIT-Sciensano/RiSpy/releases/download/v1.0.0/rispy_data.tar.gz")

# All files that MUST exist after extraction
REQUIRED_FILES = {
    "GCF_034140825.1_ASM3414082v1_genomic.fna",
    "clusters.json",
    "hq_genotype.pkl",
    "chromosome_mapping.txt",
    "3kall_snpposition_map.tsv",
    "3kall_variety_map.tsv",
    "3K-HDRA-snp-comm-miss5pc_index_position.tsv",
    "representatives.fasta",
    "CHAINFILE.chain",
}

DEEPVARIANT_CHUNKS = [
    "deepvariant.tar.gz.part-aa",
    "deepvariant.tar.gz.part-ab",
    "deepvariant.tar.gz.part-ac",
    "deepvariant.tar.gz.part-ad",
    "deepvariant.tar.gz.part-ae",
]

DEEPVARIANT_URL_BASE = ("https://github.com/BIOIT-Sciensano/RiSpy/releases/download/v1.0.0/")

# --- DeepVariant venv settings ---
DEEPVARIANT_VERSION = "1.10.0"
DEEPVARIANT_PYTHON = "3.10"          # the Python the DeepVariant 1.10 binaries/packages are built for
DEEPVARIANT_TOOLS = ("make_examples", "call_variants", "postprocess_variants")
VENV_MARKER = ".rispy_venv_ok"       # written only after the venv is fully built and verified
SKIP_REQUIREMENTS = {"pip", "setuptools", "wheel"}
EXTRA_REQUIREMENTS = ["Pillow"]


def get_data_dir() -> Path:
    """
    Returns the path to rispy_cli/data inside the installed package.
    """
    return Path(__file__).resolve().parent / "data"

def missing_required_files() -> bool:
    """
    Returns True if ANY required file is missing from rispy_cli/data.
    """
    data_dir = get_data_dir()
    if not data_dir.exists():
        return True

    existing = {p.name for p in data_dir.glob("*")}
    return not REQUIRED_FILES.issubset(existing)

def download_and_extract():
    """
    Streams the data tarball straight from the URL and merges its contents
    into rispy_cli/data.
    """
    data_dir = get_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading and extracting data from {DATA_URL} into {data_dir}")
    with urllib.request.urlopen(DATA_URL) as response, \
            tarfile.open(fileobj=response, mode="r|gz") as tar:
        for member in tar:
            # Skip directories
            if member.isdir():
                continue
            # Strip the top-level "data/" prefix
            relative_path = Path(member.name).relative_to("data")
            # Final destination inside rispy_cli/data/
            target_path = data_dir / relative_path
            # Ensure parent directories exist
            target_path.parent.mkdir(parents=True, exist_ok=True)
            # Extract file in blocks instead of reading it whole into memory
            with tar.extractfile(member) as src, open(target_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
    print("Data download and extraction completed.")

def get_deepvariant_dir() -> Path:
    """
    Returns the path to rispy_cli/deepvariant inside the installed package.
    """
    return Path(__file__).resolve().parent / "deepvariant"

def get_deepvariant_version_dir() -> Path:
    return get_deepvariant_dir() / DEEPVARIANT_VERSION

def get_deepvariant_bin_dir() -> Path:
    return get_deepvariant_version_dir() / "bazel-out" / "k8-opt" / "bin" / "deepvariant"

def get_deepvariant_venv_python() -> Path:
    return get_deepvariant_version_dir() / "venv" / "bin" / "python3"

def deepvariant_binaries_ok() -> bool:
    """
    True if all DeepVariant binaries exist and are executable.
    """
    for tool in DEEPVARIANT_TOOLS:
        p = get_deepvariant_bin_dir() / tool
        if not p.exists() or not (p.stat().st_mode & stat.S_IXUSR):
            return False
    return True

def deepvariant_venv_ok() -> bool:
    """
    True if the venv was built completely, at the current location.
    The marker stores the venv python path, so moving/reinstalling the package
    triggers a rebuild.
    """
    marker = get_deepvariant_version_dir() / VENV_MARKER
    venv_python = get_deepvariant_venv_python()
    return (marker.exists()
            and marker.read_text().strip() == str(venv_python)
            and venv_python.exists())

def missing_deepvariant() -> bool:
    """
    Returns True if DeepVariant is missing, incomplete, has non-executable
    binaries, or its venv has not been built for this server.
    """
    return not (deepvariant_binaries_ok() and deepvariant_venv_ok())


def _download_deepvariant_binaries():
    """
    Downloads the chunked tarball and extracts it into rispy_cli/deepvariant.
    """
    dv_dir = get_deepvariant_dir()
    dv_dir.mkdir(parents=True, exist_ok=True)
    # Temp dir next to the install (same filesystem): HPC /tmp is often small,
    # and the final move becomes an instant rename.
    with tempfile.TemporaryDirectory(dir=dv_dir, prefix=".download_") as tmpdir:
        tmpdir = Path(tmpdir)
        staging = tmpdir / "staging"
        print("Downloading DeepVariant chunks...")
        chunk_paths = []
        for chunk in DEEPVARIANT_CHUNKS:
            dest = tmpdir / chunk
            print(f"  - {chunk}")
            urllib.request.urlretrieve(DEEPVARIANT_URL_BASE + chunk, dest)
            chunk_paths.append(str(dest))
        print("Extracting DeepVariant...")
        cat = subprocess.Popen(["cat", *chunk_paths], stdout=subprocess.PIPE)
        with tarfile.open(fileobj=cat.stdout, mode="r|gz") as tar:
            for member in tar:
                # Skip directories and symlinks
                if member.isdir() or member.issym() or member.islnk():
                    continue
                # Strip the top-level "deepvariant/" prefix
                try:
                    relative_path = Path(member.name).relative_to("deepvariant")
                except ValueError:
                    continue
                target_path = staging / relative_path
                target_path.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(member) as src, open(target_path, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                # Make the main DeepVariant binaries executable
                if relative_path.name in DEEPVARIANT_TOOLS:
                    target_path.chmod(0o755)
        cat.stdout.close()
        if cat.wait() != 0:
            raise RuntimeError("Merging DeepVariant chunks with `cat` failed.")
        # Move the finished extraction into place
        for item in staging.iterdir():
            final = dv_dir / item.name
            if final.exists():
                shutil.rmtree(final)
            item.rename(final)

def _find_uv() -> str:
    """
    Locates the uv executable: the `uv` pip package first, then PATH.
    """
    try:
        from uv import find_uv_bin
        return find_uv_bin()
    except Exception:
        pass
    uv = shutil.which("uv")
    if uv:
        return uv
    raise RuntimeError("uv is required to build the DeepVariant environment. "
                       "Install it with: pip install uv")

def _write_requirements(version_dir: Path) -> Path:
    """
    Pins the packages of the venv shipped in the tarball into
    dv_requirements.txt (reused on later rebuilds).
    """
    req_file = version_dir / "dv_requirements.txt"
    if req_file.exists():
        return req_file
    pins = []
    for d in sorted(version_dir.glob("venv/lib/python3*/site-packages/*.dist-info")):
        name, version = d.name[:-len(".dist-info")].rsplit("-", 1)
        if name.lower() not in SKIP_REQUIREMENTS:
            pins.append(f"{name}=={version}")
    if not pins:
        raise RuntimeError(f"No packages found in the shipped venv under {version_dir}.")
    req_file.write_text("\n".join(pins) + "\n")
    return req_file

def _patch_launcher(launcher: Path, python_path: Path):
    """
    DeepVariant's Bazel launchers are zip files that hard-code the build
    machine's interpreter (PYTHON_BINARY). Point it at the venv python.
    Entries are streamed one by one, so the launcher is never loaded whole.
    """
    tmp = launcher.with_name(launcher.name + ".patching")
    with open(launcher, "rb") as fin:
        zin = zipfile.ZipFile(fin)
        preamble_len = min(i.header_offset for i in zin.infolist())
        fin.seek(0)
        preamble = fin.read(preamble_len)          # the shebang line before the zip
        main = zin.read("__main__.py").decode()
        new_main, n = re.subn(r"^PYTHON_BINARY = .*$",
                              f"PYTHON_BINARY = '{python_path}'", main, count=1, flags=re.M)
        if n != 1:
            raise RuntimeError(f"No PYTHON_BINARY line found in {launcher.name}.")
        with open(tmp, "wb") as fout:
            fout.write(preamble)
            with zipfile.ZipFile(fout, "w") as zout:
                for info in zin.infolist():
                    if info.filename == "__main__.py":
                        zout.writestr(info, new_main.encode())
                    else:
                        with zin.open(info) as src, zout.open(info, "w") as dst:
                            shutil.copyfileobj(src, dst)
    shutil.copymode(launcher, tmp)
    os.replace(tmp, launcher)

def _build_deepvariant_venv():
    """
    Replaces the shipped venv with a fresh one:
    a standalone Python 3.10 stored inside the DeepVariant folder, the same
    package versions as the shipped venv, and launchers patched to use it.
    """
    version_dir = get_deepvariant_version_dir()
    venv_dir = version_dir / "venv"
    venv_python = get_deepvariant_venv_python()
    (version_dir / VENV_MARKER).unlink(missing_ok=True)
    uv = _find_uv()
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("VIRTUAL_ENV", None)
    env["UV_PYTHON_INSTALL_DIR"] = str(version_dir / "python")      # Python lives inside DV
    env["UV_PYTHON_BIN_DIR"] = str(version_dir / "python" / "bin")  # don't touch ~/.local/bin
    env["UV_CACHE_DIR"] = str(version_dir / ".uv_cache")            # don't fill the home quota
    def run(*cmd):
        subprocess.run([uv, *cmd], check=True, env=env)
    req_file = _write_requirements(version_dir)
    print(f"Building DeepVariant environment (Python {DEEPVARIANT_PYTHON})...")
    if venv_dir.exists():
        shutil.rmtree(venv_dir)
    run("python", "install", DEEPVARIANT_PYTHON)
    run("venv", "--seed", "--python-preference", "only-managed",
        "--python", DEEPVARIANT_PYTHON, str(venv_dir))
    run("pip", "install", "--python", str(venv_python), "-r", str(req_file),
        *EXTRA_REQUIREMENTS)
    shutil.rmtree(version_dir / ".uv_cache", ignore_errors=True)
    print("Pointing DeepVariant launchers at the new environment...")
    for tool in DEEPVARIANT_TOOLS:
        _patch_launcher(get_deepvariant_bin_dir() / tool, venv_python)
    # Sanity check before declaring success: start every tool so that ALL its
    print("Checking that all DeepVariant tools start...")
    for tool in DEEPVARIANT_TOOLS:
        result = subprocess.run([str(get_deepvariant_bin_dir() / tool), "--helpshort"],
                                capture_output=True, text=True, env=env)
        output = result.stdout + result.stderr
        if "Traceback" in output:
            raise RuntimeError(f"DeepVariant `{tool}` fails to start in the new "
                               f"environment:\n{output[-2000:]}")
    (version_dir / VENV_MARKER).write_text(str(venv_python) + "\n")

def download_and_extract_deepvariant():
    """
    Installs DeepVariant: downloads the binaries only if they are missing,
    then (re)builds the venv if needed.
    """
    if not deepvariant_binaries_ok():
        _download_deepvariant_binaries()
    if not deepvariant_venv_ok():
        _build_deepvariant_venv()
    print("DeepVariant installation is done.")