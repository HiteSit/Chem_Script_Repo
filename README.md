# Contents
- [Dock_Score_AdGPU](https://github.com/HiteSit/Chem_Script_Repo/tree/master/Dock_Score_AdGPU): OLD
- [Gold_Classic_Scaffold_Docking](https://github.com/HiteSit/Chem_Script_Repo/tree/master/Gold_Classic_Scaffold_Docking): OLD
- [GOLD_Reverse_Docking](https://github.com/HiteSit/Chem_Script_Repo/tree/master/GOLD_Reverse_Docking): OLD

- [OpenMM_Tleap_PlainMD](https://github.com/HiteSit/Chem_Script_Repo/tree/master/OpenMM_Tleap_PlainMD)
- [General_Utils](https://github.com/HiteSit/Chem_Script_Repo/tree/master/General_Utils)

## OpenMM Tleap
This repository includes a tleap script to generate the topology and parameter files for OpenMM simulations.

The script initiates with a call to tleap to prepare the topology file. This process is fully automated, ensuring the topology, the simulation box size, and the ion concentration are smoothly established. The forcefield applied to the macromolecule can be customized within the tleap section of the script, as well as that for the ligand.

After topology generation, the script facilitates a classical molecular dynamics simulation using a Langevin integrator. The equilibration phase consists of both NVT (constant number of particles, volume, and temperature) and NPT (constant number of particles, pressure, and temperature) ensembles. During the NVT phase, restraints are automatically applied to the solute (excluding the ions) to stabilize the water molecules. For the NPT phase, restraints can be specified through a dictionary, allowing for selective constraints on either macromolecular atoms (e.g., the backbone) or entire residues, such as the ligand.

Finally, a standard molecular dynamics simulation is conducted at 300K.

## GeneralUtils/Fix_Protonate.py
This Python library provides a set of tools for the preparation, protonation, and standardization of small molecules. It combines the functionality of several cheminformatics libraries, including RDKit, OpenEye, and DataMol, to facilitate molecular preprocessing for computational chemistry and molecular modeling tasks.

## Features

- **Canonical SMILES Generation**: Generate canonical SMILES representations of molecules, with options for including isomeric and kekulized forms.
- **Protonation**: Generate all possible protonation states (protomers) for a given molecule.
- **Standardization**: Standardize molecular structures to ensure consistency in representation.
- **3D Structure Generation**: Generate 3D coordinates for molecules, suitable for use in molecular dynamics or docking studies.
- **Preprocessing for Searching**: Prepare molecules for database searching by standardizing and reionizing them.

## Usage
```python
from Fix_Protonate import prepare_small_mol, prepare_for_searching

smile = "CC(=O)Oc1ccccc1C(=O)O"
mol_prep_3d, smile_prep = prepare_small_mol(smile, gen_3d=True, ID="TMP", protonate=True)

```
