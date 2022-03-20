import itertools
from typing import Optional, List

from rdkit import Chem as rdChem
import numpy as np

def GetBestConformerRMS(
    mol: rdChem.Mol,
    heavyAtomsOnly: bool = False,
    confIds: Optional[List[int]] = None
) -> np.ndarray:
    if confIds is None:
        confIds = [c.GetId() for c in mol.GetConformers()]
    n_conformers = len(confIds)

    if heavyAtomsOnly:
        mol = rdChem.RemoveHs(mol)

    rms = np.zeros((n_conformers, n_conformers))
    for i, j in itertools.combinations(np.arange(n_conformers), 2):
        rmsd = rdChem.rdMolAlign.GetBestRMS(mol, mol, confIds[i], confIds[j])
        rms[i, j] = rms[j, i] = rmsd
    return rms