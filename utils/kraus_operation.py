from __future__ import annotations

import numpy as np

from pauli_string import PauliString
from pauli_operator import PauliOperator

class KrausOperation:

    def __init__(
            self,
            operators: list[PauliOperator],
    ) -> None:
        self._operators = operators
        self._adjoint_operators = [op.adjoint for op in operators]
        self._chi_matrix = None
        self._pauli_transfer_matrix = None

    def __str__(self) -> str:
        ...

    def __call__(self, operator: PauliOperator | PauliString) -> PauliOperator:
        if isinstance(operator, PauliString):
            operator = PauliOperator({operator.normalized: operator.coeff})

        res = 0
        for op, op_dag in zip(self.operators, self.adjoint_operators):
            res += op * operator * op_dag
        return res

    @property
    def operators(self) -> list[PauliOperator]:
        return self._operators

    @property
    def adjoint_operators(self) -> list[PauliOperator]:
        return self._adjoint_operators

    @property
    def chi_matrix(self) -> np.ndarray:
        ...

    @property
    def choi_matrix(self) -> np.ndarray:
        ...

    @property
    def pauli_transfer_matrix(self) -> np.ndarray:
        ...

    @property
    def is_diagonal(self) -> bool:
        return all([len(op.strings) == 1 for op in self.operators])

    def _operators_to_chi_matrix(self) -> np.ndarray:
        ...

    def _operators_to_ptm(self) -> np.ndarray:
        ...

if __name__ == "__main__":
    dm = PauliOperator(
        {PauliString.from_string("I"): 0.5,
         PauliString.from_string("Z"): -0.5}
    )

    print(dm)
    print(dm.matrix_representation)

    p = 0.1
    bitflip = KrausOperation(
        [
            PauliOperator({PauliString.from_string("I"): np.sqrt(1.0 - p)}),
            PauliOperator({PauliString.from_string("x"): np.sqrt(p)})
        ]
    )

    dm_flip = bitflip(dm)
    print(dm_flip)
    print(dm_flip.matrix_representation)

