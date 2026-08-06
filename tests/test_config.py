from __future__ import annotations

from pathlib import Path

from experimental_glycosylation_sites.config import load_config


def write_config(tmp_path: Path, target_name: str) -> Path:
    (tmp_path / "config").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sibling").mkdir(parents=True, exist_ok=True)
    (tmp_path / "sibling" / "real.csv").write_text("a\n")
    cfg = tmp_path / "config" / "test.toml"
    cfg.write_text(
        "[paths]\n"
        f'pairs_master = "../sibling/{target_name}"\n'
        'cache_dir = "cache"\n'
        "[layers]\n"
        "uniprot = true\n"
        "[policy]\n"
        'qualifying_uniprot_tiers = ["manual_experimental"]\n'
        "[api]\n"
        "delay_seconds = 0.4\n"
    )
    return cfg


def test_paths_resolve_relative_to_config_file_not_cwd(tmp_path, monkeypatch):
    cfg = write_config(tmp_path, "real.csv")
    monkeypatch.chdir(tmp_path.parent)
    config = load_config(cfg)
    assert config.paths["pairs_master"] == (tmp_path / "sibling" / "real.csv").resolve()


def test_validate_inputs_reports_missing_file(tmp_path):
    cfg = write_config(tmp_path, "absent.csv")
    errors = load_config(cfg).validate_inputs()
    assert len(errors) == 1
    assert "pairs_master" in errors[0]
    assert "absent.csv" in errors[0]


def test_validate_inputs_empty_when_all_present(tmp_path):
    cfg = write_config(tmp_path, "real.csv")
    assert load_config(cfg).validate_inputs() == []


def test_config_hash_is_stable_and_content_sensitive(tmp_path):
    first = load_config(write_config(tmp_path, "real.csv")).config_hash
    again = load_config(write_config(tmp_path / "b", "real.csv")).config_hash
    different = load_config(write_config(tmp_path / "c", "other.csv")).config_hash
    assert first == again
    assert first != different


def test_output_directories_resolve_to_the_module_root(tmp_path):
    """cache_dir and results_dir must land beside src/, not inside config/."""
    module_root = Path(__file__).resolve().parents[1]
    config = load_config(module_root / "config" / "default.toml")
    assert config.paths["cache_dir"] == module_root / "data" / "cache"
    assert config.paths["results_dir"] == module_root / "results"
