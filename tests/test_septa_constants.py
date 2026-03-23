import importlib.util
from pathlib import Path


def load_septa_module():
    repo_root = Path(__file__).resolve().parents[1]
    septa_file = repo_root / "Resonant Extraction example" / "septa.py"
    spec = importlib.util.spec_from_file_location("septa_module", septa_file)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_septa_name_groups_have_expected_sizes():
    septa = load_septa_module()

    assert len(septa.septum_namesZS) == 5
    assert len(septa.septum_namesMST) == 3
    assert len(septa.septum_namesMSE) == 5


def test_positions_vectors_are_consistent():
    septa = load_septa_module()

    assert len(septa.ZSentrances) == len(septa.ZSexits) == len(septa.blade_positionsZS)
    assert len(septa.MSTentrances) == len(septa.MSTexits) == len(septa.blade_positionsMST)
    assert len(septa.blade_positionsMSEentrance) == len(septa.blade_positionsMSEexits)


def test_critical_parameters_are_positive_and_ordered():
    septa = load_septa_module()

    assert septa.p > 0
    assert septa.Brho > 0
    assert septa.apertureZS > 0
    assert septa.apertureMST > 0
    assert septa.apertureMSE > 0

    assert all(start < end for start, end in zip(septa.ZSentrances, septa.ZSexits))
    assert all(start < end for start, end in zip(septa.MSTentrances, septa.MSTexits))
