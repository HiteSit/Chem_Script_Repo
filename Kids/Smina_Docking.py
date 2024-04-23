import os
import subprocess
from tempfile import gettempdir
import argparse

import datamol as dm
from rdkit import Chem
from rdkit.Chem import AllChem

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from pymol import cmd
from openbabel import pybel

def pybel_converter(input, input_format, output, output_format):
    mols = list(pybel.readfile(input_format, input))
    out = pybel.Outputfile(output_format, output, overwrite=True)

    for mol in mols:
        # Fix kekulization
        mol.OBMol.PerceiveBondOrders()
        mol.OBMol.AssignSpinMultiplicity(True)
        mol.OBMol.FindRingAtomsAndBonds()

        # Add Charges
        charges = mol.calccharges(model="gasteiger")

        out.write(mol)
    out.close()

class Moloc_Vina:
    def __init__(self):
        pass

    @staticmethod
    def set_vina_env():
        os.environ['PATH'] = '/home/hitesit/Software/ADFRsuite_x86_64Linux_1.0/bin:' + os.environ.get('PATH', '')
        os.environ["MGLROOT"] = "/home/hitesit/Software/mgltools_x86_64Linux2_1.5.7"

    @staticmethod
    def shape_box(selection='sele', extending=6.0, software='vina'):

        ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)

        minX = minX - float(extending)
        minY = minY - float(extending)
        minZ = minZ - float(extending)
        maxX = maxX + float(extending)
        maxY = maxY + float(extending)
        maxZ = maxZ + float(extending)

        SizeX = maxX - minX
        SizeY = maxY - minY
        SizeZ = maxZ - minZ
        CenterX = (maxX + minX) / 2
        CenterY = (maxY + minY) / 2
        CenterZ = (maxZ + minZ) / 2

        cmd.delete('all')

        if software == 'vina':
            return {'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ}, {'size_x': SizeX, 'size_y': SizeY,
                                                                                     'size_z': SizeZ}
        elif software == 'ledock':
            return {'minX': minX, 'maxX': maxX}, {'minY': minY, 'maxY': maxY}, {'minZ': minZ, 'maxZ': maxZ}
        elif software == 'both':
            return ({'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ},
                    {'size_x': SizeX, 'size_y': SizeY, 'size_z': SizeZ}), (
            {'minX': minX, 'maxX': maxX}, {'minY': minY, 'maxY': maxY}, {'minZ': minZ, 'maxZ': maxZ})

    def _getbox(self, protein, ligand):
        # Load the protein and crystal ligand
        cmd.reinitialize()
        cmd.load(protein, object="protein_XXX")
        cmd.load(ligand, object="ligand_XXX")

        # Compute the binding site
        center, size = self.shape_box(selection='ligand_XXX', extending=6.0, software='vina')
        center = list(center.values())
        size = list(size.values())
        return center, size

    @staticmethod
    def _prepare_protein(protein_to_prep):
        protein_basename = os.path.basename(protein_to_prep).split(".")[0]
        protein_prepared = f"{protein_basename}_prepared.pdb"

        # PDBFixer
        fixer = PDBFixer(filename=protein_to_prep)
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        PDBFile.writeFile(fixer.topology, fixer.positions, open(protein_prepared, 'w'))
        return protein_prepared

    @staticmethod
    def choose_protein(protein_prep):
        return os.path.abspath(protein_prep)

    @staticmethod
    def _prepare_ligand(sdf_file, ligand_name):

        # If the ligand is SDF
        if os.path.exists(sdf_file):
            # Load and add hydrogen
            mol = dm.read_sdf(sdf_file, sanitize=False)[0]

            ############################################
            # Back and forward
            tmp_smile = dm.to_smiles(mol, kekulize=True)
            mol = dm.to_mol(tmp_smile)
            ############################################

            mol_h = dm.add_hs(mol, add_coords=True)
            mol_h.SetProp("_Name", ligand_name)

            AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())

            ligand_prepared_path = os.path.join(gettempdir(), f"{ligand_name}_prepared.sdf")
            dm.to_sdf(mol_h, ligand_prepared_path)
            return ligand_prepared_path

        # If the ligand is SMILES
        else:
            mol = dm.to_mol(sdf_file)
            mol_h = Chem.AddHs(mol, addCoords=True)
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
            # Set the name
            mol_h.SetProp("_Name", ligand_name)
            # Save the ligand
            ligand_prepared_path = os.path.join(gettempdir(), f"{ligand_name}_prepared.sdf")
            dm.to_sdf(mol_h, ligand_prepared_path)
            return ligand_prepared_path

    @staticmethod
    def runner(ligand_prepared, ligand_name, protein_prepared, crystal_ligand, mode):

        docked_ligand = f"{ligand_name}_Docked.sdf"

        if mode == "docking":
            runner_command = f"smina -r {protein_prepared} -l {ligand_prepared} --autobox_ligand {crystal_ligand} -o {docked_ligand} --log {ligand_name}.log"
            subprocess.run(runner_command, shell=True, check=True, text=False)
        elif mode == "rescoring":
            runner_command = f"smina -r {protein_prepared} -l {ligand_prepared} --autobox_ligand {crystal_ligand} -o {docked_ligand} --minimize"
            subprocess.run(runner_command, shell=True, check=True, text=False)

        return docked_ligand

    def run_docking(self, ligand_to_dock, ligand_name, protein_to_dock, crystal_ligand, mode):
        protein_prepared_path = self._prepare_protein(protein_to_dock)
        ligand_prepared_path = self._prepare_ligand(ligand_to_dock, ligand_name)

        docked_ligand_path = self.runner(ligand_prepared_path,
                                         ligand_name=ligand_name,
                                         protein_prepared=protein_prepared_path,
                                         crystal_ligand=crystal_ligand,
                                         mode=mode
                                         )
        return docked_ligand_path

@cmd.extend
def pymol_docking(ligand_selection, ligand_name, protein_selection, mode):

    tmp_save_path_ligand = os.path.join(gettempdir(), f"{ligand_name}.sdf")
    tmp_save_path_protein = os.path.join(gettempdir(), "myprot.pdb")

    # Select and save the ligand to dock
    cmd.select(ligand_name, ligand_selection)
    cmd.save(tmp_save_path_ligand, ligand_name)

    pybel_converter(tmp_save_path_ligand, "sdf", tmp_save_path_ligand, "sdf")

    # Select and save the protein to dock
    cmd.select("protein", protein_selection)
    cmd.save(tmp_save_path_protein, "protein")

    print("Running Docking - Do not TOUCH the GUI")
    docker = Moloc_Vina()
    docked_ligand_path = docker.run_docking(ligand_to_dock=tmp_save_path_ligand,
                                           ligand_name=ligand_name,
                                           protein_to_dock=tmp_save_path_protein,
                                           crystal_ligand='crystal.sdf',
                                           mode=mode
                                           )
    cmd.load(docked_ligand_path)
    os.remove("./myprot_prepared.pdb")