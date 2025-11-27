import numpy as np

from pauli_string import PauliString

class PauliOperator:

    _matrix_representation: np.ndarray = None

    def __init__(
            self,
            string_coeff_pairs: dict[PauliString: complex]
    ) -> init:
        self._string_coeff_pairs = string_coeff_pairs

    @property
    def string_coeff_pairs(self) -> dict[PauliString: complex]:
        return self._string_coeff_pairs

    @property
    def strings(self) -> list[PauliString]:
        return list(self._string_coeff_pairs.keys())

    @property
    def coeffs(self) -> list[PauliString]:
        return list(self._string_coeff_pairs.values())

    @property
    def matrix_representation(self):
        res = 0
        for P, C in self.string_coeff_pairs:
            res = C * P.matrix_representation + res
        return res

class KrausOperation:

    _operators: list[PauliOperator] | None = None
    _chi_matrix: np.ndarray | None = None
    _pauli_transfer_matrix: np.ndarray | None = None

    def __init__(
            self,
            *args,
            operators: list[PauliOperator] | None = None,
            chi_matrix: np.ndarray | None = None,
            pauli_transfer_matrix: np.ndarray | None = None
    ) -> None:
        if operators is None and chi_matrix is None and pauli_transfer_matrix is None:
            raise ValueError(f"Please initialize using one of: operators, chi_matrix or pauli_transfer_matrix")