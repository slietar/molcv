from typing import Sequence
import numpy as np

from . import _molcv


def compute_cv(residue_atom_counts: Sequence[int] | np.ndarray, atom_positions: Sequence[Sequence[float]] | np.ndarray, *, cutoffs: Sequence[float] | np.ndarray):
  atom_positions_ = np.asarray(atom_positions, dtype=np.float32)
  cutoffs_ = np.ascontiguousarray(cutoffs, dtype=np.float32)
  residue_atom_counts_ = np.ascontiguousarray(residue_atom_counts, dtype=np.uint32)

  if atom_positions_.ndim != 2:
    raise ValueError('Atom positions must be a 2D array')

  if not (3 <= atom_positions_.shape[1] <= 4) :
    raise ValueError('Atom positions must have 3 or 4 columns')

  if atom_positions_.shape[1] == 3:
    atom_positions_ = np.c_[atom_positions_, np.zeros(atom_positions_.shape[0], dtype=np.float32)]

  if residue_atom_counts_.ndim != 1:
    raise ValueError('Atom counts must be a 1D array')

  if cutoffs_.ndim != 1:
    raise ValueError('Cutoffs must be a 1D array')

  # print(bytes(memoryview(atoms.astype(np.uint8))))

  atom_positions_ = np.ascontiguousarray(atom_positions_)

  return _molcv.compute_cv(residue_atom_counts_, atom_positions_, cutoffs_)
