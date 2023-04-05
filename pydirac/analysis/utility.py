import re

from pydirac.core.periodic import Element

__all__ = ["get_keyword", "get_orbital_info", "get_energy", "get_symbol_and_charge"]


def get_keyword(element: str, task_type: str, orbital_info: str, delimer: str = "@") -> str:
    """
    Constructs a keyword for a given element, task type, and orbital information.

    Parameters:
    -----------
    element : str
        The element symbol for which the keyword is to be constructed.
    task_type : str
        The type of task associated with the keyword.
    orbital_info : str
        The orbital information to be included in the keyword.
    delimer : str, optional
        The delimiter to be used to join the different components of the keyword.

    Returns:
    --------
    keyword : str
        The constructed keyword.

    Example:
    --------
    >>> get_keyword('C', 'energy', '(core 2)[vir 3]')
    'C@energy@(core 2)[vir 3]'
    """
    return delimer.join([Element(element).symbol, task_type, orbital_info])


def get_orbital_info(nocc: int, nvir: int) -> str:
    """
    Constructs the orbital information string for a given number of occupied and virtual orbitals.

    Parameters:
    -----------
    nocc : int
        The number of occupied orbitals.
    nvir : int
        The number of virtual orbitals.

    Returns:
    --------
    orbital_info : str
        The constructed orbital information string.

    Example:
    --------
    >>> get_orbital_info(2, 3)
    '(core 2)[vir 3]'
    """
    return f"(core {nocc})[vir {nvir}]"


def get_symbol_and_charge(filename="atom.mol"):
    """
    Get atom symbols and charges from a `.mol` file.

    Parameters:
    -----------
    filename : str, optional
        The name of the `.mol` file to read.

    Returns:
    --------
    atoms : list of tuple of str and float
        A list of tuples containing the symbol and charge of each atom in the file.

    Example:
    --------
    >>> get_symbol_and_charge('atom.mol')
    [('C', 6.0), ('H', 1.0), ('H', 1.0), ('H', 1.0)]
    """

    with open(filename) as f:
        context = f.read()

    # Regular expression to match the symbol and charge of each atom in the file
    pattern = r"^\s+(\d+)\.\s+(\d+\.?\d?)\s+"
    re_obj = re.compile(pattern)
    atoms = re.findall(pattern, context, re.MULTILINE)

    # Convert charge from string to float and return a list of (symbol, charge) tuples
    atoms = [(re.sub(r"[0-9]+", "", symbol), float(charge)) for symbol, charge in atoms]

    return atoms


def get_energy(filename, method="CCSD(T)"):
    """
    Get the energy from a quantum chemical calculation output file.

    Given a quantum chemical calculation output file specified by `filename` and
    the type of energy to retrieve given by `method`, the options are "CCSD(T)",
    "CCSD", "MP2", or "SCF".

    Parameters:
    -----------
    filename : str
        The name of the quantum chemical calculation output file.
    method : str, optional
        The type of energy to retrieve. Default is "CCSD(T)".

    Returns:
    --------
    energy : float
        The energy value.

    Raises:
    -------
    RuntimeError :
        If the energy for the specified method is not found in the file.

    Example:
    --------
    >>> get_energy('calc.out', method='CCSD(T)')
    -245.186792800
    """

    with open(filename) as f:
        context = f.read()

    # Regular expression to match the energy for the specified method
    if "(" and ")" in method:
        # Escape parentheses
        method = method.replace("(", r"\(").replace(")", r"\)")
    pattern = r"^@.*{method} energy[\s:]*(-?[\d\.]+)".format(**{"method": method})
    energy = re.findall(pattern, context, re.MULTILINE)

    if energy:
        return float(energy[-1])
    else:
        raise RuntimeError(f"Did not find energy for {method} in {filename}")
