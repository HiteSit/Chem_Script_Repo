# Gold Classic Scaffold Docking
This repository contains scripts for two distinct docking approaches using GOLD:
- Virtual Screening Docking using an SDF file (created with DataWarrior).
- Scaffold Docking using an SDF file (created with DataWarrior), biased through a predefined scaffold.

# OpenMM Tleap
This repository includes a tleap script to generate the topology and parameter files for OpenMM simulations.

The script initiates with a call to tleap to prepare the topology file. This process is fully automated, ensuring the topology, the simulation box size, and the ion concentration are smoothly established. The forcefield applied to the macromolecule can be customized within the tleap section of the script, as well as that for the ligand.

After topology generation, the script facilitates a classical molecular dynamics simulation using a Langevin integrator. The equilibration phase consists of both NVT (constant number of particles, volume, and temperature) and NPT (constant number of particles, pressure, and temperature) ensembles. During the NVT phase, restraints are automatically applied to the solute (excluding the ions) to stabilize the water molecules. For the NPT phase, restraints can be specified through a dictionary, allowing for selective constraints on either macromolecular atoms (e.g., the backbone) or entire residues, such as the ligand.

Finally, a standard molecular dynamics simulation is conducted at 300K.