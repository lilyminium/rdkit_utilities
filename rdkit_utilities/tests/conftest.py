import pytest
from rdkit_utilities import rdmolops


@pytest.fixture
def propylparaben():
    from rdkit_utilities import rdmolfiles
    rdmol = rdmolfiles.MolFromSmiles(
        ("[H:1]-[O:2]-[c:3]1:[c:4](-[H:5]):[c:6](-[H:7])"
         ":[c:8](-[C:13](=[O:14])-[O:15]-[C:16](-[H:17])"
         "(-[H:18])-[C:19](-[H:20])(-[H:21])-[C:22](-[H:23])"
         "(-[H:24])-[H:25]):[c:9](-[H:10]):[c:11]:1-[H:12]"),
        removeHs=False
    )
    rdmol = rdmolops.OrderByMapNumber(rdmol, clearAtomMapNumbers=True)
    return rdmol
