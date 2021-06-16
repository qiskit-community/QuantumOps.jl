using ElectronicStructure
using QuantumOps
using ZChop

## using PyCall will trigger loading pyscf-specific code

using PyCall

using ElectronicStructure: one_electron_integrals,
   MolecularData, PySCF

##  Define geometries: specification of atoms and their positions
##  for various molecules.

geoms = (
    Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414))),

    Geometry(Atom(:Li, (0., 0., 0.)), Atom(:H, (0., 0., 1.4))),

    Geometry(Atom(:O, (0., 0., 0.)), Atom(:H, (0.757, 0.586, 0.)),
             Atom(:H, (-0.757, 0.586, 0.)))
)

## Choose one of the geometries
geom = geoms[1]

#basis = "sto-3g"
basis = "631g"

## Construct specification of electronic structure problem
mol_spec = MolecularSpec(geometry=geom, basis=basis)

## Do calculations and populate MolecularData with results
mol_pyscf = MolecularData(PySCF, mol_spec)

zchop!(mol_pyscf.one_body_integrals)
zchop!(mol_pyscf.two_body_integrals)

## Create interaction operator from one and two body integrals and constant.
## This does just a bit of manipulation of mol_data; converting space orbitals
## into space-and-spin orbitals.
## This is the same as the operator by the same name in OpenFermion.
## !!! NOTE Qiskit order is: index_order=:physicist, block_spin=true
iop = InteractionOperator(mol_pyscf) # ; index_order=:default, block_spin=true)

fop = FermiSum(iop)

nothing
