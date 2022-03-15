"""
rdkit_utilities.py
RDKit utilities

"""

from rdkit import Chem as rdChem
import numpy as np

def MolFromSmiles(
    smiles: str,
    allowCXSMILES: bool = True,
    maxIterations: int = 0,
    parseName: bool = False,
    removeHs: bool = True,
    sanitize: bool = True,
    strictCXSMILES: bool = True,
    useLegacyStereo: bool = True,
) -> rdChem.Mol:
    """
    Create a molecule from SMILES.

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
    return rdChem.MolFromSmiles(smiles, smiles_parser)



def OrderByMapNumber(
    mol: rdChem.Mol,
    clearAtomMapNumbers: bool = True,
) -> rdChem.Mol:
    """
    Reorder RDKit molecule by atom map number

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        RDKit molecule
    clearAtomMapNumbers: bool
        Whether to set atom map numbers in the output
        molecule to 0

    Returns
    -------
    rdkit.Chem.Mol
        Output molecule. This is a **copy** of the input molecule.

    Example
    -------

    ::

        >>> rdmol = Chem.MolFromSmiles("[C:3][C:2][O:1]")
        >>> rdmol.GetAtomWithIdx(0).GetSymbol()
        'C'
        >>> reordered = OrderByMapNumber(rdmol)
        >>> reordered.GetAtomWithIdx(0).GetSymbol()
        'O'
    """
    map_numbers = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    order = list(map(int, np.argsort(map_numbers)))
    reordered = rdChem.RenumberAtoms(mol, order)
    if clearAtomMapNumbers:
        for atom in reordered.GetAtoms():
            atom.SetAtomMapNum(0)
    return reordered
