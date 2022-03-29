import re
import pytest

import numpy as np
from numpy.testing import assert_allclose

from rdkit import Chem as rdChem
from rdkit_utilities import rdchem, rdmolfiles
from rdkit_utilities.tests.conftest import propylparaben


@pytest.mark.parametrize(
    "smarts, indices",
    [
        ("[C](=[O])-[O-]", [(1, 2)]),
        (
            "[N:1](-[H:2])(-[H:3])-[C:4](-[O:5]-[H:6])(-[O:7]-[H:8])-[O:9]-[H:10]",
            [(1, 2), (4, 6, 8), (5, 7, 9)],
        ),
        (
            "[N-2:1]-[C:2](-[O:3]-[H:4])(-[O:5]-[H:6])-[O:7]-[H:8]",
            [(2, 4, 6), (3, 5, 7)],
        ),
    ],
)
def test_GetSymmetricAtomIndices(smarts, indices):
    mol = rdmolfiles.MolFromSmarts(
        smarts, orderByMapNumber=True, clearAtomMapNumbers=True
    )
    rdChem.SanitizeMol(mol)
    assert rdchem.GetSymmetricAtomIndices(mol) == indices


def test_AddConformerWithCoordinates(propylparaben):
    assert propylparaben.GetNumConformers() == 0
    coordinates = np.ones((25, 3))
    rdchem.AddConformerWithCoordinates(propylparaben, coordinates)
    assert propylparaben.GetNumConformers() == 1
    rdcoords = propylparaben.GetConformer(0).GetPositions()
    assert_allclose(np.array(rdcoords), coordinates)


def test_AddConformerWithCoordinates_error(propylparaben):
    assert propylparaben.GetNumConformers() == 0
    coordinates = np.ones((24, 3))

    err = "Shape of coordinates must be (25, 3). " "Given array with shape (24, 3)"

    with pytest.raises(ValueError, match=re.escape(err)):
        rdchem.AddConformerWithCoordinates(propylparaben, coordinates)
    assert propylparaben.GetNumConformers() == 0


@pytest.mark.parametrize(
    "smarts, indices, mapping",
    [
        ("[H:1]-[O:2]", [(0, 1)], [{1: 0, 2: 1}]),
        ("[H:4]-[O:2]-[c]", [(1, 0)], [{2: 1, 4: 0}]),
        (
            "[C:2](-[H:4])(-[H:6])-[H]",
            [
                (21, 22, 23),
                (21, 23, 22),
                (21, 22, 24),
                (21, 24, 22),
                (21, 23, 24),
                (21, 24, 23),
            ],
            [
                {2: 21, 4: 22, 6: 23},
                {2: 21, 4: 23, 6: 22},
                {2: 21, 4: 22, 6: 24},
                {2: 21, 4: 24, 6: 22},
                {2: 21, 4: 23, 6: 24},
                {2: 21, 4: 24, 6: 23},
            ],
        ),
    ],
)
def test_GetTaggedSubstructMatches(propylparaben, smarts, indices, mapping):
    query = rdChem.MolFromSmarts(smarts)
    as_list = rdchem.GetTaggedSubstructMatches(propylparaben, query)
    assert set(as_list) == set(indices)

    as_dict = rdchem.GetTaggedSubstructMatches(
        propylparaben, query, mappingAsTaggedDict=True
    )
    assert len(as_dict) == len(mapping)
    for mapping_ in mapping:
        assert mapping_ in as_dict


def test_SetPropsFromDict():
    properties = dict(
        boolean_prop1=True,
        boolean_prop2=False,
        integer_prop1=0,
        integer_prop2=1,
        integer_prop3=-10,
        np_int_prop=np.int64(4),
        float_prop=3.0,
        np_float_prop=np.float64(12),
        str_prop="asdf",
        empty_str_prop="",
    )

    atom = rdChem.Atom("H")
    rdchem.SetPropsFromDict(atom, properties)
    atom_properties = atom.GetPropsAsDict()
    assert properties == atom_properties
    for key in properties:
        assert np.issubdtype(type(properties[key]), type(atom_properties[key]))

    mol = rdChem.MolFromSmiles("CC")
    rdchem.SetPropsFromDict(mol, properties)
    mol_properties = mol.GetPropsAsDict()
    assert properties == mol_properties
    for key in properties:
        assert np.issubdtype(type(properties[key]), type(mol_properties[key]))
