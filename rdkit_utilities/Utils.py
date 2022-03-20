import functools
from typing import Optional, List

import numpy as np


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
            from .rdmolops import OrderByMapNumber
            mol = OrderByMapNumber(mol, clearAtomMapNumbers=clearAtomMapNumbers)
        elif clearAtomMapNumbers:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
        return mol
    return wrapper


def compute_atom_distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    dist_sq = np.einsum('ijk,ilk->ijl', coordinates, coordinates)
    diag = np.einsum("ijj->ij", dist_sq)
    a, b = diag.shape
    dist_sq += dist_sq - diag.reshape((a, 1, b)) - diag.reshape((a, b, 1))
    diag[:] = -0.0
    return np.sqrt(-dist_sq)


def get_maximally_diverse_indices(
    distance_matrix: np.ndarray,
    distance_threshold: float = 0.05,
    n_indices: Optional[int] = None,
) -> List[int]:
    n_distances = len(distance_matrix)
    if distance_matrix.shape != (n_distances, n_distances):
        raise ValueError("`distance_matrix` should be square distance matrix")

    if n_indices is None:
        n_indices = n_distances
    n_indices = min(n_indices, n_distances)

    selected_indices = [0]
    for i in range(n_indices - 1):
        selected_rms = distance_matrix[selected_indices]
        any_too_close = np.any(selected_rms < distance_threshold, axis=0)
        if np.all(any_too_close):
            break

        rmsdist = np.where(any_too_close, -np.inf, selected_rms.sum(axis=0))
        selected_indices.append(rmsdist.argmax())

    return selected_indices
