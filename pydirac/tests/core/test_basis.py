import os
import tempfile

from pydirac.core.basis import get_custom_basis_from_ele


def test_basis():
    filename = "B_c-dyall.acv4z.mol"
    with tempfile.TemporaryDirectory() as tmpdirname:
        filepath = os.path.join(tmpdirname, filename)
        get_custom_basis_from_ele("B", "acv4z", filepath)
