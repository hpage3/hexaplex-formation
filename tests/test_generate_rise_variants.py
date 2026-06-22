import importlib.util
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "generate_rise_variants.py"
SPEC = importlib.util.spec_from_file_location("generate_rise_variants", SCRIPT_PATH)
generate_rise_variants = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = generate_rise_variants
SPEC.loader.exec_module(generate_rise_variants)


def test_model_id_for_rise_values():
    assert generate_rise_variants.model_id_for(30, 3.40) == "twist30_rise3p40"
    assert generate_rise_variants.model_id_for(30, 3.38) == "twist30_rise3p38"


def test_parse_rise_values():
    assert generate_rise_variants.parse_rise_values(["3.40", "3.38"]) == [3.40, 3.38]


def test_rise_token_roundtrip_style():
    assert generate_rise_variants.rise_token(3.38) == "3p38"
