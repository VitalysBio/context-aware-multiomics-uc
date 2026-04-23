from pathlib import Path


def get_project_root() -> Path:
    """
    Return the project root directory assuming this file lives in:
    src/multiomics_uc/paths.py
    """
    return Path(__file__).resolve().parents[2]


def get_path_from_root(*parts: str) -> Path:
    """
    Build a path relative to the project root.
    """
    return get_project_root().joinpath(*parts)
