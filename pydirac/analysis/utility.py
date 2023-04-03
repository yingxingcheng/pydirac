import re


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

    with open(filename, "r") as f:
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

    with open(filename, "r") as f:
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
        raise RuntimeError("Did not find energy for {0} in {1}".format(method, filename))
