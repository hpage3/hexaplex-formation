import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


audit = _load_script_module("audit_hxc590_s1_pnab_yaml_inputs", "scripts/audit_hxc590_s1_pnab_yaml_inputs.py")


def write_yaml_fixture(path: Path, backbone: str = "backbone.pdb", base: str = "base.pdb") -> None:
    path.write_text(
        "\n".join(
            [
                "Base 1:",
                "  code: O",
                f"  file_path: {base}",
                "Backbone:",
                f"  file_path: {backbone}",
                "HelicalParameters:",
                "  h_rise: [3.4, 3.4, 1]",
                "  rise: [3.4, 3.4, 1]",
                "  h_twist: [24, 24, 1]",
                "  twist: [24, 24, 1]",
                "RuntimeParameters:",
                "  is_hexad: true",
                "  build_strand: [true, true, true, true, true, true]",
                "  strand_orientation: [true, false, true, false, true, false]",
                "  strand: OOOOOOOOOO",
                "  num_candidates: 1",
                "  num_steps: 10000000",
                "  ff_type: MMFF94",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_yaml_audit_extracts_hexad_fields(tmp_path):
    yaml_path = tmp_path / "Hexad_Antiparallel.yaml"
    (tmp_path / "backbone.pdb").write_text("PDB\n", encoding="utf-8")
    (tmp_path / "base.pdb").write_text("PDB\n", encoding="utf-8")
    write_yaml_fixture(yaml_path)

    row = audit.audit_yaml(yaml_path, tmp_path)

    assert row["h_rise"] == "3.4;3.4;1"
    assert row["h_twist"] == "24;24;1"
    assert row["is_hexad"] == "True"
    assert row["strand_orientation"] == "True;False;True;False;True;False"
    assert row["backbone_exists"] == "yes"
    assert row["base_files_exist"] == "yes"


def test_missing_file_references_are_reported(tmp_path):
    yaml_path = tmp_path / "Hexad.yaml"
    write_yaml_fixture(yaml_path, backbone="missing_backbone.pdb", base="missing_base.pdb")

    row = audit.audit_yaml(yaml_path, tmp_path)

    assert row["backbone_exists"] == "no"
    assert row["base_files_exist"] == "no"


def test_generic_yaml_is_not_hxc590_without_provenance(tmp_path):
    yaml_path = tmp_path / "Hexad.yaml"
    write_yaml_fixture(yaml_path)

    row = audit.audit_yaml(yaml_path, tmp_path)

    assert row["classification"] == "generic pNAB example YAML"
    assert "no HXC590 provenance marker" in row["provenance_note"]


def test_smoke_failure_is_captured_without_crashing(monkeypatch, tmp_path):
    monkeypatch.setattr(audit, "pnab_import_status", lambda: (False, "python cannot import pnab"))
    rows = [{"path": str(tmp_path / "Hexad.yaml")}]

    smoke = audit.smoke_rows(rows)

    assert smoke[0]["status"] == "not_run_external_tool_not_available"
    assert "cannot import pnab" in smoke[0]["message"]
    assert "not installed or integrated" in smoke[0]["message"]


def test_report_contains_conservative_provenance_language(monkeypatch, tmp_path):
    monkeypatch.setattr(audit, "pnab_import_status", lambda: (False, "python cannot import pnab"))
    report = tmp_path / "report.md"
    yaml_rows = [
        {
            "path": str(tmp_path / "Hexad.yaml"),
            "classification": "generic pNAB example YAML",
            "backbone_exists": "yes",
            "base_files_exist": "yes",
            "h_rise": "3.4;3.4;1",
            "h_twist": "24;24;1",
            "is_hexad": "True",
            "strand_orientation": "True;False",
            "num_steps": "10000000",
        }
    ]
    smoke = audit.smoke_rows(yaml_rows)

    audit.write_report(report, tmp_path, yaml_rows, smoke)

    text = report.read_text(encoding="utf-8")
    assert "pNAB exists as a separate local utility repository" in text
    assert "does not vendor, wrap, install, or integrate pNAB" in text
    assert "no pNAB-generated twist variants should be added to this branch" in text
