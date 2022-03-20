from typing import Optional, Dict

from rdkit import Chem as rdChem


from rdkit_utilities.utils import reorder_constructed_molecule

__all__ = [
    "MolFromSmiles",
    "MolFromSmarts",
]


@reorder_constructed_molecule
def MolFromSmiles(
    smiles: str,
    allowCXSMILES: bool = True,
    maxIterations: int = 0,
    parseName: bool = False,
    removeHs: bool = True,
    sanitize: bool = True,
    strictCXSMILES: bool = True,
    useLegacyStereo: bool = True,
    orderByMapNumber: bool = False,
    clearAtomMapNumbers: bool = False,
) -> Optional[rdChem.Mol]:
    """
    Create a molecule from a SMILES string.

    Example
    -------

    ::

        >>> rdmol = MolFromSmiles("[C:3][C:2][H:1]")
        >>> rdmol.GetNumAtoms()
        2
        >>> rdmol = MolFromSmiles("[C:3][C:2][H:1]", removeHs=False)
        >>> rdmol.GetNumAtoms()
        3
    """
    smiles_parser = rdChem.rdmolfiles.SmilesParserParams()
    smiles_parser.allowCXSMILES = allowCXSMILES
    smiles_parser.maxIterations = maxIterations
    smiles_parser.parseName = parseName
    smiles_parser.removeHs = removeHs
    smiles_parser.sanitize = sanitize
    smiles_parser.strictCXSMILES = strictCXSMILES
    smiles_parser.useLegacyStereo = useLegacyStereo
    mol = rdChem.MolFromSmiles(smiles, smiles_parser)
    return mol


@reorder_constructed_molecule
def MolFromSmarts(
    smarts: str,
    mergeHs: bool = False,
    replacements: Dict[str, str] = {},
    orderByMapNumber: bool = False,
    clearAtomMapNumbers: bool = False,
) -> Optional[rdChem.Mol]:
    """
    Create a molecule from a SMARTS string.
    """
    mol = rdChem.MolFromSmarts(smarts,
                               mergeHs=mergeHs,
                               replacements=replacements)
    return mol
