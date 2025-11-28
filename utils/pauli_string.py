from __future__ import annotations

import numpy as np

I = np.array([[1,0],[0,1]])
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])

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

        self._stabilizer_dict = {}
        self._matrix_representation = None
    def __str__(self) -> str:
        string_coeff = np.real_if_close(self.coeff)
        if self.coeff == 1:
            res = ""
        else:
            res = f"{string_coeff} * "
        for i in range(self._qubit_width):
            if i in self.xs:
                res += f"X({i})"
            elif i in self.ys:
                res += f"Y({i})"
            elif i in self.zs:
                res += f"Z({i})"
        return res

    def __mul__(self, other: PauliString) -> PauliString:
        new_coeff = self.coeff * other.coeff
        num_qubits = max([self.qubit_width, other.qubit_width])
        stab_1 = self.stabilizer(num_qubits)
        stab_2 = other.stabilizer(num_qubits)
        new_stab = stab_1 ^ stab_2

        zs_1, xs_1 = np.split(stab_1, 2)
        zs_2, xs_2 = np.split(stab_2, 2)
        pow_vector = (zs_1 * xs_2 - xs_1 * zs_2)*(1 - 2*zs_1*xs_1)*(1 - 2*zs_2*xs_2)
        int_coeff = (1j)**np.sum(pow_vector)
        return PauliString.from_stabilizer(new_stab, coeff=new_coeff*int_coeff)

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
        return self._ys

    @property
    def zs(self) -> tuple[int]:
        return self._zs

    @property
    def matrix_representation(self) -> np.ndarray:
        if self._matrix_representation is None:
            self._matrix_representation = np.array([[1]], dtype = np.complex128)
            for i in range(self.qubit_width):
                if i in self.xs:
                    self._matrix_representation = np.kron(X,self._matrix_representation)
                elif i in self.ys:
                    self._matrix_representation = np.kron(Y,self._matrix_representation)
                elif i in self.zs:
                    self._matrix_representation = np.kron(Z,self._matrix_representation)
                else:
                    self._matrix_representation = np.kron(I,self._matrix_representation)
            self._matrix_representation *= self.coeff
        return self._matrix_representation

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

        return PauliString(tuple(xs),tuple(zs),tuple(ys),coeff=coeff, qubit_width=len(pauli_string))

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

        res = PauliString(xs,zs,coeff=coeff, qubit_width=num_qubits)
        res._stabilizer_dict = {-1:stabilizer}
        return res

    def stabilizer(self, num_qubits: int = -1, swap: bool = False) -> np.ndarray:
        if -1 < num_qubits < self.qubit_width:
            raise ValueError(f"num_qubits({num_qubits}) must be greater than qubit_width({self.qubit_width})")
        if num_qubits not in self._stabilizer_dict:
            if num_qubits == -1:
                self._stabilizer_dict[num_qubits] = self._make_stabilizer()
            else:
                self._stabilizer_dict[num_qubits] = self._make_stabilizer(num_qubits)

        if swap:
            return np.concatenate(
                np.split(self._stabilizer_dict[num_qubits], 2)[::-1]
            )
        return self._stabilizer_dict[num_qubits]

    def _make_stabilizer(self, num_qubits: int | None = None) -> np.ndarray:
        """
        Generate the stabilizer representation of the pauli string
        """
        if num_qubits is None:
            num_qubits = self._qubit_width

        zs_stab = np.zeros(num_qubits, dtype = int)
        xs_stab = np.zeros(num_qubits, dtype = int)

        for i in self.xs:
            xs_stab[i] ^= 1
        for i in self.zs:
            zs_stab[i] ^= 1
        for i in self.ys:
            xs_stab[i] ^= 1
            zs_stab[i] ^= 1

        return np.concatenate((zs_stab, xs_stab))

    def symplectic_product(self, other: PauliString) -> int:
        """
        Compute the symplectic product of two pauli strings
        """
        num_qubits = max(self.qubit_width, other.qubit_width)
        stab_1 = self.stabilizer(num_qubits)
        stab_2 = other.stabilizer(num_qubits, swap=True)
        return sum(stab_1 & stab_2) % 2

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
