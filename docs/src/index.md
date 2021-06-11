# QuantumOps.jl

Documentation for QuantumOps.jl

This package implements types representing single-particle Pauli and fermionic operators and
their composition into terms and sums.

Types implemented include single-qubit(fermion) Pauli(fermionic) opearators,
mulit-qubit Pauli(fermionic) operators, and linear combinations of multi-particle operators

```@autodocs
Modules = [QuantumOps.AbstractOps]
Order   = [:function, :type]
```

```@docs
OpTerm
```

```@docs
OpTerm(x)
```
