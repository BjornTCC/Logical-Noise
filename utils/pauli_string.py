from __future__ import annotations

import numpy as np

class PauliString:

    def __init__(
        self,
        xs: tuple[int] = (),
        zs: tuple[int] = (),
        ys: tuple[int] | None = None,
        coeff: int | float | None = None,
        qubit_width: int | None = None
    ) -> None:

        if ys is None:
            ys = tuple([i for i in set(xs).intersection(set(zs))])
            xs = tuple([i for i in xs if i not in ys])
            zs = tuple([i for i in zs if i not in ys])

        self.xs = xs
        self.zs = zs
        self.ys = ys

        if coeff is not None:
            self._coeff = coeff
        else:
            self._coeff: float = 1

        if qubit_width is None:
            self._qubit_width = max(list(xs) + list(ys) + list(zs)) + 1
        else:
            self._qubit_width = qubit_width
    @property
    def qubit_width(self) -> int:
        return self._qubit_width

    @classmethod
    def from_string(cls, string: str) -> PauliString:
        print(test)

    @classmethod
    def from_stabilizer(cls, stabilizer: list[int], coeff: int | float | None = None) -> PauliString:
        if len(stabilizer) % 2 != 0:
            raise ValueError(f"Length of stabilizer must be even, got {len(stabilizer)}")

        num_qubits = len(stabilizer) // 2

        zs, xs = (),()



    def __str__(self) -> str:
        if self._coeff == 1:
            res = ""
        else:
            res = f"{self._coeff} * "

        for i in range(self._qubit_width):
            if i in self.xs:
                res += f"X({i})"
            elif i in self.ys:
                res += f"Y({i})"
            elif i in self.zs:
                res += f"Z({i})"

        return res

    def to_stabilizer(self, num_qubits: int | None = None, as_ints: bool = False, reverse: bool = False) -> list[int]:
        if num_qubits is None:
            num_qubits = self._qubit_width
        if reverse:
            res =  [(i in self.xs or i in self.ys) for i in range(num_qubits)] + [(i in self.zs or i in self.ys) for i in range(num_qubits)]
        else:
            res = [(i in self.zs or i in self.ys) for i in range(num_qubits)] + [(i in self.xs or i in self.ys) for i in
                                                                                 range(num_qubits)]
        if as_ints:
            return np.array([int(i) for i in res])
        return res

    def symplectic_product(self, other: PauliString) -> int:
        num_qubits = max(self.qubit_width, other.qubit_width)
        stab_1 = self.to_stabilizer(num_qubits, as_ints=True)
        stab_2 = other.to_stabilizer(num_qubits, as_ints=True, reverse=True)
        return sum(stab_1 & stab_2)

    def commutative_sign(self, other: PauliString) -> int:
        return (-1)**self.symplectic_product(other)

if __name__ == "__main__":
    print("test")
    tuples = (
        (0,2),
        (1,2),
    )

    pstring = PauliString(*tuples)
    other = PauliString((0,))
    print(pstring)
    print(other)
    print(pstring.commutative_sign(other))