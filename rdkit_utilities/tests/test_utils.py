import pytest

from rdkit_utilities import Chem, Utils


@pytest.fixture
def ordered_molecule():
    rdmol = Chem.MolFromSmiles(
        ("[H:1]-[O:2]-[c:3]1:[c:4](-[H:5]):[c:6](-[H:7])"
         ":[c:8](-[C:13](=[O:14])-[O:15]-[C:16](-[H:17])"
         "(-[H:18])-[C:19](-[H:20])(-[H:21])-[C:22](-[H:23])"
         "(-[H:24])-[H:25]):[c:9](-[H:10]):[c:11]:1-[H:12]"),
        removeHs=False
    )
    rdmol = Chem.OrderByMapNumber(rdmol, clearAtomMapNumbers=True)
    return rdmol


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
    ordered_molecule,
    core_indices,
    include_central_atoms,
    n_neighbors,
    expected_indices
):
    indices = Utils.GetAtomNeighborIndices(
        ordered_molecule,
        centralAtomIndices=core_indices,
        includeCentralAtoms=include_central_atoms,
        numAtomNeighbors=n_neighbors
    )
    assert indices == expected_indices
