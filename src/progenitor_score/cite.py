"""Attribution text shipped with the tool (also in ATTRIBUTION.md, CITATION.cff and model.json)."""
ATLAS = ("Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, "
         "Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of "
         "aberrant AML differentiation. Blood Cancer Discov 2025;6:307-24. doi:10.1158/2643-3230.BCD-24-0342")
TOOL = "Conor Finlay Lab, Trinity College Dublin. progenitor-score (code MIT; model tables CC BY-NC 4.0)."
SHORT = "progenitor-score: model trained on the BoneMarrowMap atlas, Zeng et al. Blood Cancer Discov 2025;6:307-24 (doi:10.1158/2643-3230.BCD-24-0342); model tables CC BY-NC 4.0."

def citation() -> str:
    """The two references to cite when you use this tool."""
    return f"Please cite:\n  1. {ATLAS}\n  2. {TOOL}\nModel tables: CC BY-NC 4.0 (LICENSE-MODEL.md)."
