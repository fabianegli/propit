# © Fabian Egli, FGCZ, ETHZ
# GPLv3+

import csv
import logging
import os
from pathlib import Path

import pandas as pd

from .constant import CRUX_COMET_PIN_COLUMN_DTYPES

logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
logger.addHandler(console_handler)


def read_pin_tab(pin_file):
    """Parse a Percolator input (.pin) file into a pandas DataFrame.

    Thie pin file format allows for variable lenght lines because if a peptide is
    assigned to mupltiple proteins, these identifiers are also separated by tabs,
    the same separator normally used for column separation.

    NB: This parser will combine the multiple protein identifiers per PSM by
    concatenation with a semicolon (;).

    For more information about the FlashLFQ generic input format plese consult
    https://github.com/smith-chem-wisc/FlashLFQ/wiki/Identification-Input-Formats

    Args:
        pin_file: Path to the .pin file.

    Returns:
        pandas.DataFrame: A DataFrame containing the parsed data from the .pin file.
    """

    rows = []
    with open(pin_file, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")

        headers = next(reader)
        n_static_columns = len(headers) - 1

        for row in reader:
            fixed = row[:n_static_columns]
            protein_ids = row[n_static_columns:]
            rows.append(fixed + [";".join(protein_ids)])

    return pd.DataFrame(
        rows,
        columns=headers,
    ).astype({k: v for k, v in CRUX_COMET_PIN_COLUMN_DTYPES.items() if k in headers})


def read_psms_tab(psms_file):
    return pd.read_table(psms_file, low_memory=False)


def read_comet_pin_combined(comet_output_dir: Path):
    comet_output_dir = Path(comet_output_dir)
    pin_files = sorted(comet_output_dir.glob("*.pin"))
    if not pin_files:
        raise ValueError(f"'.pin' files are required. Found none in {comet_output_dir}")

    dfs = []

    for pin_file in pin_files:
        logger.info(f"Reading: {pin_file}")
        dfs.append(read_pin_tab(pin_file=pin_file))

    df = pd.concat(dfs, ignore_index=True)

    return df


def read_comet_txt_combined(comet_output_dir: Path, crux=False):
    comet_output_dir = Path(comet_output_dir)
    if not comet_output_dir.is_dir():
        raise ValueError(f"{comet_output_dir=} must be a directory.")

    decoy_files = sorted(comet_output_dir.glob("*.decoy.txt"))

    if not decoy_files:
        raise ValueError("Decoys are required.")

    dfs = []

    for decoy_file in decoy_files:
        stem = decoy_file.name.removesuffix(".decoy.txt")
        target_file = comet_output_dir / (stem + ".txt")
        decoys = pd.read_table(decoy_file, low_memory=False, skiprows=1).assign(Label=-1)
        targets = pd.read_table(target_file, low_memory=False, skiprows=1).assign(Label=1)
        logger.info(f"Reading {decoy_file} ...")
        dfs.append(decoys)
        dfs.append(targets)

    df = pd.concat(dfs, ignore_index=True).assign(
        file=stem,
        SpecId=lambda x: (
            x.file.str.cat(x.scan.astype("str"), sep="_")
            .str.cat(x.charge.astype("str"), sep="_")
            .str.cat(x.num.astype("str"), sep="_")
        ),
    )

    return df


def read_percolator_target(target_file):
    return pd.read_table(target_file, low_memory=False).assign(PSMId=lambda x: x.PSMId.map(os.path.basename))
