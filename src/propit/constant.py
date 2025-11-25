# © Fabian Egli, FGCZ, ETHZ
# GPLv3+

FLASHLFQ_GENERIC_INPUT_COLUMNS = [
    "File Name",
    "Scan Retention Time",
    "Precursor Charge",
    "Base Sequence",
    "Full Sequence",
    "Peptide Monoisotopic Mass",
    "Protein Accession",
]

PROTEOBENCH_GENERIC_UPLOAD_COLUMNS = [
    "Sequence",
    "Proteins",
    "Charge",
    "Modified sequence",
    "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01",
    "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02",
    "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03",
    "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01",
    "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02",
    "LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03",
]

CRUX_COMET_PIN_COLUMN_DTYPES = {
    "Label": "int",
    "ScanNr": "int",
    "rt": "float",
    "ExpMass": "float",
    "CalcMass": "float",
    "deltLCn": "float",
    "deltCn": "float",
    "XCorr": "float",
    "TailorScore": "float",
    "byIonsMatched": "int",
    "byIonsTotal": "int",
    "byIonsFraction": "float",
    "byIonsRepeatMatch": "int",
    "PepLen": "int",
    "Charge1": "int",
    "Charge2": "int",
    "Charge3": "int",
    "Charge4": "int",
    "Charge5": "int",
    "enzN": "int",
    "enzC": "int",
    "enzInt": "int",
    "lnNumDSP": "float",
    "dM": "float",
    "absdM": "float",
    "score": "float",
    "q-value": "float",
    "posterior_error_prob": "float",
}
