import pytest

from pydiscamb import get_table

@pytest.mark.parametrize(
        [
            "table_str",
            "expected_entries",
            "expected_number_of_gaussian_parameters"
        ],
        [
            (
                "Waasmeier-Kirfel",
                211,
                5,
            ),
            (
                "electron-cctbx",
                98,
                5,
            ),
            (
                "electron-IT",
                98,
                5,
            ),
            (
                "IT92",
                213,
                4,
            ),
        ]
)
def test_get_table(table_str: str, expected_entries: int, expected_number_of_gaussian_parameters: int):
    table = get_table(table_str)
    assert len(table) == expected_entries
    for key, val in table.items():
        assert len(val.a) == expected_number_of_gaussian_parameters
        assert len(val.b) == expected_number_of_gaussian_parameters
    