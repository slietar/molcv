import sys
from typing import Iterable, Optional, overload

import numpy as np

from . import _molcv


def cli():
  _molcv.cli(sys.argv)


@overload
def compute_cv(
  residue_atom_counts: Iterable[int],
  atom_positions: Iterable[Iterable[float]],
  *,
  cutoff: float
) -> np.ndarray:
  ...

@overload
def compute_cv(
  residue_atom_counts: Iterable[int],
  atom_positions: Iterable[Iterable[float]],
  *,
  cutoffs: Iterable[float]
) -> np.ndarray:
  ...

def compute_cv(
    residue_atom_counts: Iterable[int],
    atom_positions: Iterable[Iterable[float]],
    *,
    cutoff: Optional[float] = None,
    cutoffs: Optional[Iterable[float]] = None
  ):
  """
  Calculate the circular variance of protein residues.

  Parameters
    atom_positions: The positions of the N atoms of the protein, as an array with shape (N, 3) or (N, 4) (in which case the fourth coordinate is unused).
    cutoff/cutoffs: The M cutoff distance(s) for the calculation.
    residue_atom_counts: The number of atoms in each of the P residues.

  Returns
    The circular variance of each residue as an array with shape (M, P), or (P) if `cutoff` was used.
  """

  match (cutoff, cutoffs):
    case (None, None):
      raise ValueError('Either "cutoff" or "cutoffs" must be provided')
    case (None, _):
      cutoffs_ = np.ascontiguousarray(cutoffs, dtype=np.float32)
    case (_, None):
      cutoffs_ = np.array([cutoff], dtype=np.float32)
    case (_, _):
      raise ValueError('Only one of "cutoff" or "cutoffs" can be provided')

  atom_positions_ = np.asarray(atom_positions, dtype=np.float32)
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

  atom_positions_ = np.ascontiguousarray(atom_positions_)

  result = _molcv.compute_cv(residue_atom_counts_, atom_positions_, cutoffs_)

  return result[0, :] if cutoff is not None else result


__all__ = [
  'cli',
  'compute_cv'
]
