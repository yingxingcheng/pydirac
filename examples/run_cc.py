import tempfile

from pydirac import AtomicCCSet, AtomicDHFSet, Molecule

molecule = Molecule(["Mc"], [[0.0, 0.0, 0.0]])
dhf = AtomicDHFSet(
    molecule=molecule,
    is_ff=True,
    ff_mode="Q",
    hamiltonian_mode="2C",
    is_spinfree=True,
    user_mol_settings={"basis_type": "dyall.cv3z"},
)
with tempfile.TemporaryDirectory() as tmp_dir:
    dhf.write_input(output_dir=tmp_dir)
    print(" Inp file ".center(80, "#"))
    print(dhf.inp)
    print(" Mol file ".center(80, "#"))
    print(dhf.mol)

cc = AtomicCCSet(
    molecule=molecule,
    no_T=True,
    e_min=-30,
    e_max=30,
    nelec=(2, 1),
    nel_f1=(2, 1, 0, 0),
    is_spinfree=True,
    hamiltonian_mode="4C",
    is_ff=True,
    ff_mode="D",
)

# Use a temporary directory to store the output files
with tempfile.TemporaryDirectory() as tmp_dir:
    cc.write_input(output_dir=tmp_dir)
    print(" Inp file ".center(80, "#"))
    print(cc.inp)
    print(" Mol file ".center(80, "#"))
    print(cc.mol)
