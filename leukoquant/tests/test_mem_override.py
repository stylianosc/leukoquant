"""LEUKOQUANT_TRACULA_MEM_GB: the TRACULA memory request is 8 GB unless a valid override is set."""
from leukoquant.utils import mem_override as mo

GB = 1024


def test_default_is_8_gb_when_unset_or_empty():
    assert mo.tracula_mem_mb({}) == 8 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": ""}) == 8 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "  "}) == 8 * GB


def test_valid_override_sets_whole_gb():
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "16"}) == 16 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": " 12 "}) == 12 * GB
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "4"}) == 4 * GB       # lower bound
    assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": "64"}) == 64 * GB     # upper bound


def test_invalid_override_is_ignored_with_a_warning(capsys):
    for bad in ("3", "65", "-8", "8.5", "lots", "16GB"):
        assert mo.tracula_mem_mb({"LEUKOQUANT_TRACULA_MEM_GB": bad}) == 8 * GB, bad
    assert "ignoring LEUKOQUANT_TRACULA_MEM_GB" in capsys.readouterr().out
