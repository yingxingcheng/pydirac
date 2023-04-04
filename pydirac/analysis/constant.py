from pydirac.core.periodic import Element

__all__ = ["get_keyword", "get_orbital_info"]


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
