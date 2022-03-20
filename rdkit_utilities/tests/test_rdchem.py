import pytest

from rdkit import Chem as rdChem
from rdkit_utilities import rdchem, rdmolfiles


@pytest.mark.parametrize("smarts, indices", [
    ("[C](=[O])-[O-]", [(1, 2)]),
    (
        "[N:1](-[H:2])(-[H:3])-[C:4](-[O:5]-[H:6])(-[O:7]-[H:8])-[O:9]-[H:10]",
        [(1, 2), (4, 6, 8), (5, 7, 9)]
    ),
    (
        "[N-2:1]-[C:2](-[O:3]-[H:4])(-[O:5]-[H:6])-[O:7]-[H:8]",
        [(2, 4, 6), (3, 5, 7)]
    ),
])
def test_GetSymmetricAtomIndices(smarts, indices):
    mol = rdmolfiles.MolFromSmarts(smarts, orderByMapNumber=True,
                                   clearAtomMapNumbers=True)
    rdChem.SanitizeMol(mol)
    assert rdchem.GetSymmetricAtomIndices(mol) == indices
