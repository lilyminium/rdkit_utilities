import pytest
from rdkit_utilities import Chem

@pytest.mark.parametrize("smiles, n_heavy_atoms, n_all_atoms", [
    ("[C:3][C:2][H:1]", 2, 3),
])
def test_MolFromSmiles(smiles, n_heavy_atoms, n_all_atoms):
    heavy_mol = Chem.MolFromSmiles(smiles)
    assert heavy_mol.GetNumAtoms() == n_heavy_atoms
    all_mol = Chem.MolFromSmiles(smiles, removeHs=False)
    assert all_mol.GetNumAtoms() == n_all_atoms

def test_OrderByMapNumber():
    mol = Chem.MolFromSmiles("[C:3][C:2][O:1]")
    assert mol.GetAtomWithIdx(0).GetSymbol() == "C"

    reordered = Chem.OrderByMapNumber(mol, clearAtomMapNumbers=False)
    first = reordered.GetAtomWithIdx(0)
    assert first.GetSymbol() == "O"
    assert first.GetAtomMapNum() == 1

    reordered = Chem.OrderByMapNumber(mol, clearAtomMapNumbers=True)
    first = reordered.GetAtomWithIdx(0)
    assert first.GetSymbol() == "O"
    assert first.GetAtomMapNum() == 0
    