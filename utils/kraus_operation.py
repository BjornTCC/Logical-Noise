from __future__ import annotations

import numpy as np



class KrausOperation:

    def __init__(
            self,
            *args,
            operators: list[PauliOperator],
    ) -> None:
        self._operators = operators
        self._chi_matrix = None
        self._pauli_transfer_matrix = None

    def __call__(self, operator: PauliOperator | PauliString) -> PauliOperator:
        ...

    @property
    def operators(self) -> list[PauliOperator]:
        ...

    @property
    def chi_matrix(self) -> np.ndarray:
        ...

    @property
    def pauli_transfer_matrix(self) -> np.ndarray:
        ...

    def _operators_to_chi_matrix(self) -> np.ndarray:
        ...

    def _operators_to_ptm(self) -> np.ndarray:
        ...
