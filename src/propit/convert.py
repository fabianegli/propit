# © Fabian Egli, FGCZ, ETHZ
# GPLv3+


import logging
import os
from pathlib import Path

import pandas as pd
import peptacular

from .constant import FLASHLFQ_GENERIC_INPUT_COLUMNS, PROTEOBENCH_GENERIC_UPLOAD_COLUMNS
from .read import read_comet_pin_combined, read_percolator_target, read_pin_tab, read_psms_tab

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
logger.addHandler(console_handler)


def comet_pin_to_percolator_pin(comet_output_dir: Path, pin_out_file: Path = Path("../unambiguous-pin.pin")):
    pins = read_comet_pin_combined(comet_output_dir)
    lape = pins[["Label", "Peptide"]]

    target_peptides = set(lape[lape["Label"] == 1].Peptide)
    decoy_peptides = set(lape[lape["Label"] == -1].Peptide)
    ambiguous_peptides = target_peptides.intersection(decoy_peptides)

    percolator_pin = pins[~pins.Peptide.isin(ambiguous_peptides)].drop(columns=["deltCn"])

    percolator_pin.to_csv(pin_out_file, sep="\t", index=False)
    return percolator_pin


def comet_perc_generic_flashlfqinput(
    comet_out_dir,
    percolator_output_file,
    generic_flashlfq_file,
    max_q_value=0.01,
):
    percolator_output_file = Path(percolator_output_file).absolute()
    comet_out_dir = Path(comet_out_dir).absolute()
    generic_flashlfq_file = Path(generic_flashlfq_file).absolute()

    print("Reading txt files ...")
    combined = read_comet_txt_combined(comet_out_dir)
    print("Done reading txt files.")

    print("Generate mappings ...")
    combined["file_scan"] = combined.file.str.cat(combined.scan.astype("str"), sep="_")

    map_filescan_rt_sec = make_mapping(combined, "file_scan", "retention_time_sec")
    # specid_to_rt_map = make_mapping(combined, _from="SpecId", _to="retention_time_sec")
    # specid_to_file_map = make_mapping(combined, _from="SpecId", _to="file")
    # specid_to_charge_map = make_mapping(combined, _from="SpecId", _to="charge")
    # specid_to_calcmass_map = make_mapping(combined, _from="SpecId", _to="calc_neutral_mass")
    peptide_to_calcmass_map = make_mapping(combined, _from="modified_peptide", _to="calc_neutral_mass")
    modified_to_base_peptide_map = make_mapping(combined, _from="modified_peptide", _to="plain_peptide")

    print("... finished the mappings.")

    targets = read_percolator_target(percolator_output_file)

    gen_flfq = (
        targets.query(f"`q-value` <= {max_q_value}")
        .assign(
            **{
                "PSMId": lambda x: x["PSMId"].map(os.path.basename),
                "file_scan": lambda x: x.PSMId.str.rsplit("_", n=2, expand=True).iloc[:, 0],
                "File Name": lambda x: x.PSMId.str.rsplit("_", n=3, expand=True).iloc[:, 0] + ".mzML",
                "Scan Retention Time": lambda x: x.file_scan.map(lambda x: map_filescan_rt_sec[x] / 60),
                "Precursor Charge": lambda x: x.PSMId.str.rsplit("_", n=2, expand=True).iloc[:, 1],
                "Peptide Monoisotopic Mass": lambda x: x.peptide.map(lambda x: peptide_to_calcmass_map[x]),
                "Base Sequence": lambda x: x.peptide.map(lambda x: modified_to_base_peptide_map[x]),
                "Full Sequence": lambda x: x.peptide.str[2:-2],
                "Protein Accession": lambda x: x.proteinIds,
            }
        )
        .filter(items=FLASHLFQ_GENERIC_INPUT_COLUMNS)
    )

    print("Done preparing pin data.")

    print(f"Writing generic FlashLFQ file: {generic_flashlfq_file.absolute()}")
    gen_flfq.to_csv(generic_flashlfq_file, sep="\t", index=False)

    print(f"{gen_flfq.columns=}")

    return gen_flfq


def read_comet_txt_combined(comet_output_dir: Path):
    comet_output_dir = Path(comet_output_dir)

    decoy_files = sorted(Path(comet_output_dir).glob("*.decoy.txt"))

    if not decoy_files:
        raise ValueError("Decoys are required.")

    dfs = []

    for decoy_file in decoy_files:
        print(decoy_file)
        stem = decoy_file.name.removesuffix(".decoy.txt")
        target_file = comet_output_dir / (stem + ".txt")
        decoys = pd.read_table(decoy_file, low_memory=False, skiprows=1).assign(Label=-1, file=stem)
        print(target_file)
        targets = pd.read_table(target_file, low_memory=False, skiprows=1).assign(Label=1, file=stem)
        dfs.append(decoys)
        dfs.append(targets)

    df = pd.concat(dfs, ignore_index=True).assign(
        SpecId=lambda x: (
            x.file.str.cat(x.scan.astype("str"), sep="_")
            .str.cat(x.charge.astype("str"), sep="_")
            .str.cat(x.num.astype("str"), sep="_")
        )
    )

    return df


def charge_n_to_precursor_charge(df: pd.DataFrame) -> pd.Series:
    col_charge = [(c, int(c.split("Charge")[1])) for c in df.filter(regex=r"^Charge\d+$").columns]
    return sum(df[col] * charge for col, charge in col_charge)


def specid_to_precursor_charge(df: pd.DataFrame) -> pd.Series:
    return df.SpecId.str.rsplit("_", n=2, expand=True).iloc[:, 1]


def make_mapping(df, _from, _to):
    return df[[_from, _to]].set_index(_from).to_dict()[_to]


class Percolator2FlashLFQ:
    def __init__(
        self,
        percolator_output_dir,
        max_q_value=0.01,
        remove_contaminants=True,
        cleanup=False,
        comet_perc=False,
    ):
        self.max_q_value = max_q_value
        self.remove_contaminants = remove_contaminants
        self.cleanup = cleanup
        self.comet_perc = comet_perc

        self.percolator_output_dir = Path(percolator_output_dir)
        logger.info(
            "working with percolator output files in %s",
            self.percolator_output_dir.absolute(),
        )
        self.pin_tab_file = self.percolator_output_dir / "make-pin.pin"
        self.psms_tab_file = self.percolator_output_dir / "percolator.target.psms.txt"

        self.generic_flashlfq_input = None
        self.pin_tab = None
        self.pin_tab_prepared = None
        self.psms_tab = None
        self.psms_tab_prepared = None

    def read_pin_tab(self):
        logger.info("Reading pin file ...")
        self.pin_tab = read_pin_tab(self.pin_tab_file)
        logger.info("Done reading pin file.")

    def prepare_pin_tab(self):
        self.ensure_pin_tab()

        logger.info("Preparing pin data ...")

        logger.info(f"{self.pin_tab.columns=}")

        self.pin_tab_prepared = self.pin_tab.assign(
            **{
                "Precursor Charge": charge_n_to_precursor_charge,
            }
        ).rename(
            columns={
                "CalcMass": "Peptide Monoisotopic Mass",
            }
        )
        logger.info("Done preparing pin data.")

        logger.info(f"{self.pin_tab.columns=}")

        if self.cleanup:
            del self.pin_tab

    def prepare_psms_tab(self, max_q_value, remove_contaminants):
        self.ensure_psms_tab()

        logger.info("Preparing psms data ...")

        self.max_q_value = max_q_value
        print(self.psms_tab.columns)
        self.psms_tab_prepared = (
            self.psms_tab.assign(filename=lambda x: x.PSMId.str.rsplit("_", n=2, expand=True).iloc[:, 0])
            .rename(
                columns={
                    "filename": "File Name",
                    "rt": "Scan Retention Time",
                    "proteinIds": "Protein Accession",
                }
            )
            .assign(
                **{
                    "File Name": lambda x: x["File Name"].map(os.path.basename),
                    "Full Sequence": lambda x: x.peptide.str[2:-2],
                    "Base Sequence": lambda x: x.peptide.str[2:-2],
                }
            )
            .drop(
                columns=[
                    "peptide",
                ]
            )
            .query(f"`q-value` <= {max_q_value}")
            .set_index("PSMId")
        )
        if remove_contaminants:
            self.psms_tab_prepared = self.psms_tab_prepared[
                ~self.psms_tab_prepared["Protein Accession"].str.contains("Cont_")
            ]

        n_records = self.psms_tab.shape[0]
        n_records_kept = self.psms_tab_prepared.shape[0]
        n_records_removed = n_records - n_records_kept
        logger.info(
            f"Removed {n_records_removed} psms from {n_records} leaving "
            f"{n_records_kept} with a q-value of <= {max_q_value}."
        )
        logger.info("Done preparing psms data.")

        if self.cleanup:
            del self.psms_tab

    def prepare_psms_tab_comet_percolator(self, max_q_value, remove_contaminants):
        self.ensure_psms_tab()

        logger.info("Preparing psms data ...")

        self.max_q_value = max_q_value
        logger.info(f"{self.psms_tab.columns=}")
        self.psms_tab_prepared = (
            self.psms_tab.assign(filename=lambda x: x.SpecId.str.rsplit("_", n=3, expand=True).iloc[:, 0])
            .rename(
                columns={
                    "filename": "File Name",
                    "rt": "Scan Retention Time",
                    "Proteins": "Protein Accession",
                }
            )
            .assign(
                **{
                    "File Name": lambda x: x["File Name"].map(os.path.basename),
                    "Full Sequence": lambda x: x.Peptide.str[2:-2],
                    "Base Sequence": lambda x: x.Peptide.str[2:-2],
                }
            )
            .drop(
                columns=[
                    "peptide",
                ]
            )
            .query(f"`q-value` <= {max_q_value}")
            .set_index("PSMId")
        )
        if remove_contaminants:
            self.psms_tab_prepared = self.psms_tab_prepared[
                ~self.psms_tab_prepared["Protein Accession"].str.contains("Cont_")
            ]

        n_records = self.psms_tab.shape[0]
        n_records_kept = self.psms_tab_prepared.shape[0]
        n_records_removed = n_records - n_records_kept
        logger.info(
            f"Removed {n_records_removed} psms from {n_records} leaving "
            f"{n_records_kept} with a q-value of <= {max_q_value}."
        )
        logger.info("Done preparing psms data.")

        if self.cleanup:
            del self.psms_tab

    def ensure_pin_tab(self):
        if self.pin_tab is None:
            self.read_pin_tab()

    def ensure_pin_tab_prepared(self):
        if self.pin_tab_prepared is None:
            self.prepare_pin_tab()

    def read_psms_tab(self):
        logger.info("Reading psms file ...")
        self.psms_tab = read_psms_tab(self.psms_tab_file)
        logger.info("Done reading psms file.")

    def ensure_psms_tab(self):
        if self.psms_tab is None:
            self.read_psms_tab()

    def ensure_psms_tab_prepared(self):
        if self.psms_tab_prepared is None:
            if self.comet_perc:
                prepare_psms_tab = self.prepare_psms_tab_comet_percolator
            else:
                prepare_psms_tab = self.prepare_psms_tab

        prepare_psms_tab(self.max_q_value, self.remove_contaminants)

    def copmute_generic_flashlfq_input(self):
        self.ensure_pin_tab_prepared()
        self.ensure_psms_tab_prepared()

        self.generic_flashlfq_input = self.psms_tab_prepared.join(
            self.pin_tab_prepared[self.pin_tab_prepared.Label == 1]
            .rename(columns={"SpecId": "PSMId"})
            .set_index("PSMId"),
            on="PSMId",
            how="inner",
            validate="1:1",
        )[FLASHLFQ_GENERIC_INPUT_COLUMNS]

    def ensure_generic_flashlfq_input(self):
        if self.generic_flashlfq_input is None:
            self.copmute_generic_flashlfq_input()

    def write(self, generic_flashlfq_input_file=None):
        self.ensure_generic_flashlfq_input()

        generic_flashlfq_input_file = generic_flashlfq_input_file or "generic_FlashLFQ_input-from-percolator"
        descriptive_generic_flashlfq_input_file = (
            f"{generic_flashlfq_input_file}-max_q_value{self.max_q_value}-rc{self.remove_contaminants}.tsv"
        )
        logger.info(
            "Writing FlashLFQ generic input file with shape "
            f"{self.generic_flashlfq_input.shape}: {descriptive_generic_flashlfq_input_file}"
        )
        self.generic_flashlfq_input.to_csv(descriptive_generic_flashlfq_input_file, sep="\t", index=False)


def flashlfq_peptides_to_proteobench(flashlfq_output_dir: Path):
    peptides = pd.read_table(flashlfq_output_dir / "QuantifiedPeptides.tsv", low_memory=False).assign(
        **{
            "Sequence": lambda x: x["Base Sequence"],
            "Proteins": lambda x: x["Protein Groups"],
            "Charge": 0,  # peptides are combined from charge states
            "Modified sequence": lambda x: x["Sequence"],
            "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01"
            ],
            "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02"
            ],
            "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03"
            ],
            "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01"
            ],
            "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02"
            ],
            "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03": lambda x: x[
                "Intensity_LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03"
            ],
        }
    )[PROTEOBENCH_GENERIC_UPLOAD_COLUMNS]

    # remove contaminants
    ids_before = peptides.shape[0]
    contaminants_count_before = peptides.Proteins.str.contains("Cont_").sum()
    peptides = peptides[~peptides.Proteins.str.contains("Cont_")]
    contaminants_count_after = peptides.Proteins.str.contains("Cont_").sum()
    ids_after = peptides.shape[0]

    logger.info(
        f"{contaminants_count_before} IDs from contaminants were removed from a total of {ids_before} IDs "
        f"leaving {ids_after} and {contaminants_count_after} contaminants."
    )

    proteobench_custom_upload_file = flashlfq_output_dir.parent / "ProteoBench-custom-peptidoform.tsv"

    peptides.to_csv(proteobench_custom_upload_file, sep="\t", index=False)
    logger.info(f"Generic input file for ProteoBench has been written to {proteobench_custom_upload_file.absolute()}")

    return peptides


def sage_results_2_generic_flashlfq_input(sage_output_dir, flashlfq_wd):
    sage_output_dir = Path(sage_output_dir)
    flashlfq_wd = Path(flashlfq_wd)

    if not sage_output_dir.exists():
        raise ValueError(f"{sage_output_dir=} must be an existing directory.")
    if not flashlfq_wd.exists():
        raise ValueError(f"{flashlfq_wd=} must be an existing directory.")

    results_sage = pd.read_table(sage_output_dir / "results.sage.tsv", low_memory=False)

    n_sage_results = results_sage.shape[0]
    logger.info(f"{n_sage_results} ids in results.sage.tsv")

    peptide_q_cutoff = 0.01  # noqa: F841  /  this variable is used in the pandas query

    def apply_peptide_q_cutoff(df, peptide_q_cutoff):
        logger.info(f"Applying peptide_q cut-off of {peptide_q_cutoff}")
        new = df.query("peptide_q < @peptide_q_cutoff")
        n_new = new.shape[0]
        logger.info(f"Removed {n_sage_results - n_new} leaving {n_new} records.")
        return new

    sage_flashlfq_input_preparation = (
        results_sage.pipe(apply_peptide_q_cutoff, peptide_q_cutoff=peptide_q_cutoff)
        .assign(
            **{
                "File Name": lambda x: x.filename,
                "Scan Retention Time": lambda x: x.rt,
                "Precursor Charge": lambda x: x.charge,
                "Base Sequence": lambda x: x.peptide.map(lambda x: peptacular.parse(x).sequence),
                "Full Sequence": lambda x: x.peptide,
                "Peptide Monoisotopic Mass": lambda x: x.calcmass,
                "Protein Accession": lambda x: x.proteins,
            }
        )
        .filter(items=FLASHLFQ_GENERIC_INPUT_COLUMNS)
    )

    # From https://github.com/smith-chem-wisc/FlashLFQ/wiki/Identification-Input-Formats

    # > The first line of the text file should contain column headers identifying what each column is.
    # > For search software that lists decoys and PSMs above 1% FDR, you may want to remove these
    # prior to FlashLFQ analysis.
    # > FlashLFQ will probably crash if ambiguous PSMs are passed into it (e.g., a PSM with more than
    # 2 peptides listed in one line).

    protein_ids = sage_flashlfq_input_preparation["Protein Accession"]
    decoy_sel = protein_ids.str.contains("rev_")
    logger.info(f"{decoy_sel.sum()} decoys were found.")

    full_seqs = sage_flashlfq_input_preparation["Full Sequence"]
    ambiguous_psm_sel = full_seqs.str.contains(r"\|")
    logger.info(f"{ambiguous_psm_sel.sum()} ambiguous PSMs were found.")

    if decoy_sel.sum():
        logger.info("removing decoys")
        sage_flashlfq_generic_input = sage_flashlfq_input_preparation[~decoy_sel]

        protein_ids_in = sage_flashlfq_generic_input["Protein Accession"]
        decoy_sel_in = protein_ids_in.str.contains("rev_")
        if not decoy_sel_in.sum():
            logger.info("No remaining decoys.")

    full_seqs_in = sage_flashlfq_generic_input["Full Sequence"]
    ambiguous_psm_sel_in = full_seqs_in.str.contains(r"\|")
    if not ambiguous_psm_sel_in.sum():
        logger.info("No remaining ambiguous PSMs.")

    generic_flashlfq_input_file = flashlfq_wd / "generic_FlashLFQ_input-from-sage-ion.tsv"
    sage_flashlfq_generic_input.to_csv(generic_flashlfq_input_file, sep="\t", index=False)
    logger.info(f"saved generic input for FlashLFQ to {generic_flashlfq_input_file.absolute()}")


def strip_deltacn_col(pin_file: Path, output_dir: Path | None = None):
    logger.info(f"Dropping column 'deltCn' from {pin_file.absolute()}")
    no_deltacn_dir = output_dir or pin_file.parent / "no_deltaCn"
    # TODO: also remove other columns with errors
    no_deltacn_dir.mkdir(exist_ok=True)
    stripped_pin_file = no_deltacn_dir / pin_file.name
    (
        read_pin_tab(pin_file=pin_file)
        .drop(columns=["deltCn"])
        .to_csv(
            stripped_pin_file,
            sep="\t",
            index=False,
        )
    )
    logger.info(f"Written output to new file {stripped_pin_file.absolute()}")
