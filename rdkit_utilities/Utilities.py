from typing import List, Set, Union

import numpy as np
from rdkit import Chem as rdChem


def GetAtomNeighborIndices(
    mol: rdChem.Mol,
    centralAtomIndices: Union[Set[int], List[int]] = [],
    includeCentralAtoms: bool = True,
    numAtomNeighbors: int = 0,
) -> Set[int]:
    """
    Get neighbor atom indices around specified core

    Parameters
    ----------
    mol: rdkit.Chem.Mol
    centralAtomIndices: List[int]
        Central atoms to get neighbors around
    includeCentralAtoms: bool
        Whether to include the central atoms in the output
    numAtomNeighbors: int
        Size of shell around the central core.
        -1 returns all indices in the molecule;
        0 returns no additional neighbors;
        1 returns the first shell around the central core,
        and so on.

    Returns
    -------
    Set[int]
        Set of neighbor atom indices around specified core
    """
    central_atoms = set(centralAtomIndices)
    neighbor_atoms = set(central_atoms)

    if numAtomNeighbors < 0:
        neighbor_atoms = set(range(mol.GetNumAtoms()))

    else:
        current_layer = central_atoms
        while numAtomNeighbors:
            new_layer = set()
            for index in current_layer:
                atom = mol.GetAtomWithIdx(index)
                for bond in atom.GetBonds():
                    new_layer.add(bond.GetOtherAtomIdx(index))

            new_layer -= neighbor_atoms
            neighbor_atoms |= new_layer
            current_layer = new_layer
            numAtomNeighbors -= 1

    if not includeCentralAtoms:
        neighbor_atoms -= central_atoms
    return neighbor_atoms


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
