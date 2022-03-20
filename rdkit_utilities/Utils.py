import functools


def reorder_constructed_molecule(func):
    @functools.wraps(func)
    def wrapper(
        *args,
        orderByMapNumber: bool = False,
        clearAtomMapNumbers: bool = False,
        **kwargs
    ):
        mol = func(*args, **kwargs)
        if orderByMapNumber and mol is not None:
            from .Utilities import OrderByMapNumber
            mol = OrderByMapNumber(mol, clearAtomMapNumbers=clearAtomMapNumbers)
        elif clearAtomMapNumbers:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
        return mol
    return wrapper
