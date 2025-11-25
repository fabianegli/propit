from pathlib import Path

import pandas as pd
import pytest

from propit.convert import Percolator2FlashLFQ

TEST_DATA = Path(__file__).parent / "test_data"
TEST_DATA_PERCOLATOR = TEST_DATA / "percolator-output-small"

GENERIC_FLASHLFQ_INPUT_COLUMNS = [
    "File Name",
    "Scan Retention Time",
    "Precursor Charge",
    "Base Sequence",
    "Full Sequence",
    "Peptide Monoisotopic Mass",
    "Protein Accession",
]


@pytest.fixture()
def p2f_parsed():
    p2f = Percolator2FlashLFQ(
        percolator_output_dir=TEST_DATA_PERCOLATOR, max_q_value=0.01, remove_contaminants=True, cleanup=False
    )
    p2f.copmute_generic_flashlfq_input()
    return p2f


def test_Percolator2FlashLFQ_pin_prepared_shape(p2f_parsed):
    assert p2f_parsed.pin_tab.shape == (10, 29)
    assert p2f_parsed.pin_tab_prepared.shape == (10, 30)


# @pytest.mark.parametrize("charge", [n for n in [range(1,6)]])
# this would require more careful test file comositoin to include expamples with all charge states form 2 to 5.
@pytest.mark.parametrize("charge", [n for n in [2, 3]])
def test_percolator_to_flashlfq_charge(p2f_parsed, charge):
    unique = p2f_parsed.pin_tab_prepared[p2f_parsed.pin_tab_prepared[f"Charge{charge}"] == 1][
        "Precursor Charge"
    ].unique()
    assert unique.size == 1
    assert unique[0] == charge


def test_generic_flashlfq_format(p2f_parsed):
    assert p2f_parsed.generic_flashlfq_input.columns.identical(pd.Index(GENERIC_FLASHLFQ_INPUT_COLUMNS))
