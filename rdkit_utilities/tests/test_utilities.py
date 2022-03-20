import pytest

from rdkit_utilities import Chem, Utilities


def test_OrderByMapNumber():
    mol = Chem.MolFromSmiles("[C:3][C:2][O:1]")
    assert mol.GetAtomWithIdx(0).GetSymbol() == "C"

    reordered = Utilities.OrderByMapNumber(mol, clearAtomMapNumbers=False)
    first = reordered.GetAtomWithIdx(0)
    assert first.GetSymbol() == "O"
    assert first.GetAtomMapNum() == 1

    reordered = Utilities.OrderByMapNumber(mol, clearAtomMapNumbers=True)
    first = reordered.GetAtomWithIdx(0)
    assert first.GetSymbol() == "O"
    assert first.GetAtomMapNum() == 0
    

@pytest.mark.parametrize(
    "core_indices, include_central_atoms, n_neighbors, expected_indices",
    [
        ({14}, True, -1, set(range(25))),
        ({14}, False, -1, set(range(14)) | set(range(15, 25))),
        ({14}, True, 0, {14}),
        ({14}, False, 0, set()),
        ({14}, True, 1, {12, 14, 15}),
        ({14}, False, 1, {12, 15}),
        ({14}, True, 2, {7, 13, 12, 14, 15, 16, 17, 18}),
        ({14}, False, 2, {7, 13, 12, 15, 16, 17, 18}),
        ({14}, True, 3, {5, 8, 7, 13, 12, 14, 15, 16, 17, 18, 19, 20, 21}),
        ({14}, False, 3, {5, 8, 7, 13, 12, 15, 16, 17, 18, 19, 20, 21}),
        ({2, 12, 21}, True, -1, set(range(25))),
        ({2, 12, 21}, True, 0, {2, 12, 21}),
        ({2, 12, 21}, False, 0, set()),
        ({2, 12, 21}, False, 1, {1, 3, 10, 7, 13, 14, 18, 22, 23, 24}),
        ({2, 12}, False, 2, {0, 1, 3, 4, 5, 10, 8, 11, 7, 13, 14, 5, 8, 15})
    ]
)
def test_GetAtomNeighborIndices(
    propylparaben,
    core_indices,
    include_central_atoms,
    n_neighbors,
    expected_indices
):
    indices = Utilities.GetAtomNeighborIndices(
        propylparaben,
        centralAtomIndices=core_indices,
        includeCentralAtoms=include_central_atoms,
        numAtomNeighbors=n_neighbors
    )
    assert indices == expected_indices
