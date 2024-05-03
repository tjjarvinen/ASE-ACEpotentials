# Start with
# PYTHON_JULIACALL_HANDLE_SIGNALS=yes python3

from ase_ace_test_package import ase_ace
from ase import io
# And what ever else you need from ASE

# ACEmd has the potential file and test structure.
# We need to load juliacall to get access to them.
# This is only needed for this example.
from juliacall import Main as jl
jl.seval('using ACEmd')
atoms = io.read( jl.seval('joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")') )
pot_file = jl.seval( 'joinpath(pkgdir(ACEmd), "data", "TiAl.json")' )

# Create potential
ace_pot = ase_ace.ACEcalculator(pot_file)

# and use the potential
atoms.calc = ace_pot
atoms.get_potential_energy()
atoms.get_forces()