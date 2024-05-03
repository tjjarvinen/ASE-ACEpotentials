# ASE-ACEcalculator

This python package creates ACE calculator that can used with ASE.
The calculator uses newer [ACEmd](https://github.com/ACEsuit/ACEmd.jl) implementation to perform the calculations.


## Example use

```python
from ase_ace import ase_ace
from ase import io

potential = "path to ACE potential file"

ace_pot = ase_ace.ACEpotential(potential)

atoms = io.read("atoms structure file")
atoms.calculator = ace_pot
atoms.get_forces()
```