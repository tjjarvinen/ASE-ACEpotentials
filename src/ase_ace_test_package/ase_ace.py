# You need to have PYTHON_JULIACALL_HANDLE_SIGNALS=yes environment variable set on

import numpy as np
from juliacall import Main as jl
import juliapkg
from ase import io
from ase.calculators.calculator import Calculator
from ase.constraints import voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress
#atoms = io.read( jl.seval('joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")') )
# jl.seval('using ACEmd')
# ace = ACEcalculator(pot_file)
#pot_file = jl.seval( 'joinpath(pkgdir(ACEmd), "data", "TiAl.json")' )


class ACEcalculator(Calculator):
    """
    ASE-compatible Calculator that calls ACEmd.jl for energy, forces and virial
    """
    implemented_properties = ['forces', 'energy', 'free_energy', 'stress']
    default_parameters = {}
    name = 'ACEcalculator'

    def __init__(self, potential_file, proj_dir=None):
        Calculator.__init__(self)
        if proj_dir != None:
            juliapkg.activate(proj_dir)
        jl.seval("using ACEmd.ASEhelper")
        jl.seval("using ACEmd")
        self.ace_calculator = jl.seval(f'load_ace_model("{potential_file}")') #julia.eval
        self.ace_energy = jl.seval('ase_energy')
        self.ace_forces = jl.seval('ase_forces')
        self.ace_virial = jl.seval('ase_virial')
        self.ace_energy_forces_virial = jl.seval('ase_energy_forces_virial')

    def calculate(self, atoms, properties, system_changes):
        #Calculator.calculate(self, atoms, properties, system_changes)
        atom_symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        cell = atoms.cell.cellpar()
        pbc = atoms.get_pbc()
        self.results = {}
        if 'energy' in properties and 'forces' in properties and 'stress' in properties:
            print("using combination")
            E, f, v = self.ace_energy_forces_virial(self.ace_calculator, atom_symbols, positions, cell, pbc)
            self.results['energy'] = E
            self.results['free_energy'] = E
            self.results['forces'] = np.array(f)
            voigt_stress = full_3x3_to_voigt_6_stress( v / ( -atoms.get_volume() ) )  # stress is -virial/volume
            self.results['stress'] = voigt_stress
        elif 'energy' in properties:
            E = self.ace_energy(self.ace_calculator, atom_symbols, positions, cell, pbc)
            self.results['energy'] = E
            self.results['free_energy'] = E
        elif 'forces' in properties:
            self.results['forces'] = np.array(self.ace_forces(self.ace_calculator, atom_symbols, positions, cell, pbc))
        elif 'stress' in properties:
            virial = np.array(self.ace_virial(self.ace_calculator, atom_symbols, positions, cell, pbc))
            voigt_stress = full_3x3_to_voigt_6_stress( virial / ( -atoms.get_volume() ) )  # stress is -virial/volume
            self.results['stress'] = voigt_stress
