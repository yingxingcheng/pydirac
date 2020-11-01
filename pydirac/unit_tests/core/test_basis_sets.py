from pydirac.core.basis_sets import get_custom_basis_from_ele

# get_explicit_basis_default_basis('../data/Ar_dyall.acv4z.mol', '../data/Ar_c-dyall.acv4z.mol')
# get_explicit_basis_default_basis('B', '../data/B_c-dyall.acv4z.mol')
get_custom_basis_from_ele('B', 'acv4z', '../data/B_c-dyall.acv4z.mol')