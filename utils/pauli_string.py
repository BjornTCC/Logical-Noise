from __future__ import annotations

import numpy as np

class PauliString:
    """
    Object for representing and manipulating pauli strings
    args:
        xs: tuple of qubits with X-matrix
        ys: tuple of qubits with Y-matrix
        zs: tuple of qubits with Z-matrix
        coeff: leading coefficient, if such an object is needed. Default = 1
        qubit_width: number of qubits associated to the pauli string, if None set to the maximum non-trivial qubit index
    """
    def __init__(
        self,
        xs: tuple[int] = (),
        zs: tuple[int] = (),
        ys: tuple[int] | None = None,
        coeff: int | float | complex = 1,
        qubit_width: int | None = None
    ) -> None:

        if ys is None:
            ys = tuple([i for i in set(xs).intersection(set(zs))])
            xs = tuple([i for i in xs if i not in ys])
            zs = tuple([i for i in zs if i not in ys])

        self._xs = xs
        self._zs = zs
        self._ys = ys
        self._coeff = coeff

        if qubit_width is None:
            self._qubit_width = max(list(xs) + list(ys) + list(zs)) + 1
        else:
            self._qubit_width = qubit_width

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


    @property
    def qubit_width(self) -> int:
        return self._qubit_width

    @property
    def coeff(self) -> int | float | complex | None:
        return self._coeff

    @property
    def xs(self) -> tuple[int]:
        return self._xs

    @property
    def ys(self) -> tuple[int]:
        return self._xs

    @property
    def zs(self) -> tuple[int]:
        return self._xs

    @classmethod
    def from_string(cls, string: str) -> PauliString:
        """
        Generate a PauliString from a string representation, indexing from the left
        e.g.
        -1*XIYZ -> -1 * X(0)Y(2)Z(3)
        IIZI -> Z(2)
        """
        if "*" in string:
            coeff_string, pauli_string = string.split("*")

            try:
                coeff = int(coeff_string)
            except ValueError:
                try:
                    coeff = float(coeff_string)
                except ValueError:
                    try:
                        coeff = complex(coeff_string)
                    except ValueError:
                        raise ValueError(f"Unrecognized coefficient format: {coeff_string}")
        else:
            coeff = 1
            pauli_string = string

        xs = []
        ys = []
        zs = []

        for i,pauli in enumerate(pauli_string.upper()):
            match pauli:
                case "X":
                    xs.append(i)
                case "Y":
                    ys.append(i)
                case "Z":
                    zs.append(i)
                case "I":
                    pass
                case _:
                    raise ValueError(f"Pauli strings must contain only X,Y,Z or I. Got: {pauli}")

        return PauliString(xs,zs,ys,coeff=coeff, qubit_width=len(pauli_string))

    @classmethod
    def from_stabilizer(cls, stabilizer: list[int | bool], coeff: int | float = 1) -> PauliString:
        """
        Generate from a stabilizer vector representation of the pauli string
        """
        if len(stabilizer) % 2 != 0:
            raise ValueError(f"Length of stabilizer must be even, got {len(stabilizer)}")

        num_qubits = len(stabilizer) // 2

        zs = tuple([i for i in range(num_qubits) if stabilizer[i]])
        xs = tuple([i for i in range(num_qubits) if stabilizer[num_qubits + i]])

        return PauliString(xs,zs,coeff=coeff, qubit_width=num_qubits)

    def to_stabilizer(self, num_qubits: int | None = None, as_array: bool = False, reverse: bool = False) -> list[int]:
        """
        Generate the stabilizer representation of the pauli string
        """
        if num_qubits is None:
            num_qubits = self._qubit_width
        if reverse:
            res =  [(i in self.xs or i in self.ys) for i in range(num_qubits)] + [(i in self.zs or i in self.ys) for i in range(num_qubits)]
        else:
            res = [(i in self.zs or i in self.ys) for i in range(num_qubits)] + [(i in self.xs or i in self.ys) for i in
                                                                                 range(num_qubits)]
        if as_array:
            return np.array([int(i) for i in res])
        return res

    def symplectic_product(self, other: PauliString) -> int:
        """
        Compute the symplectic product of two pauli strings
        """
        num_qubits = max(self.qubit_width, other.qubit_width)
        stab_1 = self.to_stabilizer(num_qubits, as_array=True)
        stab_2 = other.to_stabilizer(num_qubits, as_array=True, reverse=True)
        return sum(stab_1 & stab_2)

    def commutative_sign(self, other: PauliString) -> int:
        """
        Compute the sign associated with commuting two pauli strings
        """
        return (-1)**self.symplectic_product(other)

    def commutes(self, other: PauliString) -> bool:
        """
        Compute whether the pauli strings commute
        """
        return bool(self.symplectic_product(other))

if __name__ == "__main__":
    print("test")
    tuples = (
        (0,2),
        (1,2),
    )

    pstring = PauliString(*tuples)
    other = PauliString.from_string("XIYZ")
    print(pstring)
    print(other)
    print(pstring.commutative_sign(other))