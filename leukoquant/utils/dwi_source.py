"""Which TRACULA run supplies the DWI to DTI, NODDI and the metrics step.

DTI and NODDI do not depend on the parcellation, so they read one DWI per session: the one preprocessed by the TRACULA
run of the "primary" parcellation. The primary parcellation is FreeSurfer when it is in the run's list, otherwise the
first one listed. The DTI/NODDI maps already on disk were fitted on that DWI. A run that lists only another parcellation
(for example --parcellation gif) reads that parcellation's TRACULA DWI instead. The two DWI files hold the same data up to
rounding, but the metrics step amplifies the difference, so the tables of such a run are not comparable to those of a run
that lists FreeSurfer. The rule is kept, and the run states which DWI it uses.
"""
from typing import Sequence


def primary_parcellation(parcellations: Sequence[str]) -> str:
    """The parcellation whose TRACULA run supplies the DWI: freesurfer if listed, otherwise the first one."""
    return "freesurfer" if "freesurfer" in parcellations else parcellations[0]


def describe_dwi_source(parcellations: Sequence[str]) -> str:
    """One log line naming the TRACULA run the DWI is read from, with a warning when it is not the FreeSurfer one."""
    primary = primary_parcellation(parcellations)
    msg = (f"[dwi] DTI, NODDI and metrics read the DWI of tracula-{primary} "
           f"(parcellations: {', '.join(parcellations)}; freesurfer if listed, otherwise the first).")
    if primary != "freesurfer":
        msg += (" WARNING: this is not the FreeSurfer DWI. Tables from this run can differ from those of a run that lists "
                "freesurfer, because the two DWI files differ at rounding level and the metrics step amplifies that. "
                "Add freesurfer to --parcellation to get the same DWI as a full run.")
    return msg
