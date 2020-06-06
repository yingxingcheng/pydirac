from mendeleev import Element, element
from pydirac.utility.config import get_mol_by_custom_basis, \
    get_mol_by_default_basis


def get_mole_file(atom_info:Element, basis_type:str, filename_out:str,
                  basis_choice : str = 'BASIS', ) -> None:
    atom_index = atom_info.atomic_number
    atom_type = atom_info.symbol

    if basis_choice not in ['EXPLICIT', 'BASIS']:
        raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                        'for builtin basis or custom basis.')


    if basis_choice == 'EXPLICIT':
        with open('basis/{0}.dat'.format(atom_type), 'r') as f:
            basis_info = f.read()
        template = get_mol_by_custom_basis(atom_type, atom_index,
                                           basis_choice, basis_info)
    elif basis_choice == 'BASIS':
        template = get_mol_by_default_basis(atom_type, atom_index, basis_type)
    else:
        raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                        'for builtin basis or custom basis.')


    fname = filename_out or atom_type + '_' + basis_type + '.mol'
    with open(fname, 'w') as f:
        f.write(template)

