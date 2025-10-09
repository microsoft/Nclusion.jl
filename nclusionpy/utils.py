import shutil
from pathlib import Path

def clean_julia_artifacts():
    """
    Delete cached binary artifacts to force a clean re-instantiation.
    Useful on clusters when glibc/OpenSSL versions differ.
    """
    depot = Path.home() / ".julia" / "artifacts"
    if depot.exists():
        shutil.rmtree(depot)
        print(f"[nclusionpy] Cleared Julia artifacts at {depot}")
    else:
        print("[nclusionpy] No artifacts directory found.")
