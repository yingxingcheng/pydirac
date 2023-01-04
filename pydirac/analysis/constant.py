from pydirac.core.periodic import Element


def get_keyword(element, task_type, orbital_info, delimer="@"):
    return delimer.join([Element(element).symbol, task_type, orbital_info])


def get_orbital_info(nocc, nvir):
    return f"(core {nocc})[vir {nvir}]"
