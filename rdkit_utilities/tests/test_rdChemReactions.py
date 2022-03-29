from rdkit.Chem import rdChemReactions
from rdkit import Chem as rdChem
from rdkit_utilities import Chem
from rdkit_utilities.rdChemReactions import ClickReaction


def test_run_click():
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3][H]>>[C:1](=[O:2])[N:3]')
    cooh = "[H:5]-[C:1](=[O:2])-[O:3]-[H:4]"
    cnc = "[C:1](-[H:2])(-[H:3])(-[H:4])-[N:5](-[H:6])-[C:7](-[H:8])(-[H:9])(-[H:10])"
    reacts = (
        Chem.MolFromSmarts(cooh, orderByMapNumber=True, clearAtomMapNumbers=True),
        Chem.MolFromSmarts(cnc, orderByMapNumber=True, clearAtomMapNumbers=True)
    )
    products = ClickReaction(rxn, reacts)
    assert len(products) == 1

    product = products[0]
    expected = rdChem.AddHs(rdChem.MolFromSmiles("CN(C)C=O"))
    assert product.GetSubstructMatch(expected)
    assert expected.GetSubstructMatch(product)
