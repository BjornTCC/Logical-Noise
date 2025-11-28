from __future__ import annotations

import numpy as np
from copy import deepcopy

from pauli_string import PauliString

class PauliOperator:

    def __init__(
            self,
            string_coeff_pairs: dict[PauliString: complex]
    ) -> None:
        assert all([p.coeff == 1 for p in string_coeff_pairs.keys()]), f"all coefficients of keys must be 1"
        self.num_qubits = list(string_coeff_pairs.keys())[0].qubit_width
        assert all([p.qubit_width == self.num_qubits for p in string_coeff_pairs.keys()]), f"num qubits of pauli strings not all equal"
        self._string_coeff_pairs = string_coeff_pairs
        self._matrix_representation: np.ndarray = None

    def __str__(self) -> str:
        res = ""
        for p,c in self.string_coeff_pairs.items():
            if p.is_identity:
                res = f"{c}" + res
            elif isinstance(c,complex) or (c > 0):
                res += f" + {c} * {p}"
            else:
                res += f" {c} * {p}"
        return res

    def __add__(self, other: PauliOperator | PauliString | complex) -> PauliOperator:
        if isinstance(other, PauliOperator):
            new_string_coeff_pairs = deepcopy(self.string_coeff_pairs)
            for p,c in other.string_coeff_pairs.items():
                if p in new_string_coeff_pairs.keys():
                    new_string_coeff_pairs[p] += c
                else:
                    new_string_coeff_pairs[p] = c
        elif isinstance(other, PauliString):
            new_string_coeff_pairs = deepcopy(self.string_coeff_pairs)
            if other in new_string_coeff_pairs.keys():
                new_string_coeff_pairs[other] += 1.0
            else:
                new_string_coeff_pairs[other] = 1.0
        else:
            identity_string = PauliString(qubit_width=self.num_qubits)
            new_string_coeff_pairs = deepcopy(self.string_coeff_pairs)
            if identity_string in new_string_coeff_pairs.keys():
                new_string_coeff_pairs[identity_string] += other
            else:
                new_string_coeff_pairs[identity_string] == other

        return PauliOperator(new_string_coeff_pairs)

    def __radd__(self, other) -> PauliOperator:
        return self.__add__(other)

    def __mul__(self, other: PauliOperator | PauliString | complex) -> PauliOperator:
        if isinstance(other, PauliOperator):
            new_scp = {}
            for self_string, self_coeff in self.string_coeff_pairs.items():
                for other_string, other_coeff in other.string_coeff_pairs.items():
                    new_string = self_string * other_string
                    new_coeff = self_coeff * other_coeff * new_string.coeff
                    new_string = new_string.normalized

                    if new_string in new_scp.keys():
                        new_scp[new_string] += new_coeff
                    else:
                        new_scp[new_string] = new_coeff
            return PauliOperator(new_scp)

        if isinstance(other, PauliString):
            new_scp = {}
            for P, c in self.string_coeff_pairs.items():
                new_P = P * other
                new_C = c * new_P.coeff
                new_P = new_P.normalized

                if new_P in new_scp.keys():
                    new_scp[new_P] += new_C
                else:
                    new_scp[new_P] = new_C
            return PauliOperator(new_scp)

        return PauliOperator({P: c * other for P, c in self.string_coeff_pairs.items()})

    def __rmul__(self, other:  PauliString | complex) -> PauliOperator:
        if isinstance(other, PauliString):
            new_scp = {}
            for P, c in self.string_coeff_pairs.items():
                new_P =other * P
                new_C = c * new_P.coeff
                new_P = new_P.normalized
                if new_P in new_scp.keys():
                    new_scp[new_P] += new_C
                else:
                    new_scp[new_P] = new_C
            return PauliOperator(new_scp)

        return PauliOperator({P: c * other for P, c in self.string_coeff_pairs.items()})

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
        for P, C in self.string_coeff_pairs.items():
            res += C * P.matrix_representation
        return res

    @property
    def adjoint(self) -> PauliOperator:
        return PauliOperator({P: np.conjugate(C) for P, C in self.string_coeff_pairs.items()})
