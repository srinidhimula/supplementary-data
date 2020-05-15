## Rotation angle of NO2-MIL-53

Calculating the rotation angle between the benzene ring and reference plane in NO2-MIL-53 for CP2K trajectory file (.xyz)

This code uses CP2K trajectory file (.xyz), cell parameters  and number of atoms in unit cell of system being studied (in this case NO2-MIL-53) as input to calculate the rotation angle of the nitro-functionalized linker in an AIMD simulation.

- System requirements
  - Operating systems: Windows 7 or later, macOS, and Linux
  - Python 3.X version installed. 
It has been used and tested on Windows 10 and macOS Mojave.

- Instructions
  -[Rotationangle_singlecell.py](./Rotationangle_singlecell.py) can be used for a single unit cell MD simulation
  -[Rotationangle_supercell.py](./Rotationangle_supercell.py) can be used for 211 supercell MD simulation.

A demo ipython notebook file [Demo rotation angle single cell.ipynb](./DemoRotationAngleSingleCell.ipynb) is provided, which returns the rotation angle in an array for 100 timesteps from input MD trajectory file in input_demo directory (singlecell_100.xyz). These arrays can be saved into .dat files for further plotting and analysis. Expected output .dat files are also placed in (input_demo/expected_output folder)
