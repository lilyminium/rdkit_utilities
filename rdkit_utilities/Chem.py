"""
rdkit_utilities.py
RDKit utilities

"""
from typing import Optional, Dict
from typing_extensions import Literal

from rdkit import Chem as rdChem

from ._Chem import *
from ._io import ALL_RDKIT_PARSERS, molecule_from_input
from rdkit_utilities.utils import reorder_constructed_molecule


@reorder_constructed_molecule
def MolFromInput(
    molInput: str,
    *args,
    inputFormat: Optional[Literal[(*ALL_RDKIT_PARSERS,)]] = None,  # type: ignore
    orderByMapNumber: bool = False,
    clearAtomMapNumbers: bool = False,
    **kwargs,
) -> Optional[rdChem.Mol]:
    return molecule_from_input(
        molInput,
        *args,
        mol_format=inputFormat,
        **kwargs
    )
