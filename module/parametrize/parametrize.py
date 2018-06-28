from simtk.openmm import app
import openforcefield.utils as utils
from openforcefield.typing.engines.smirnoff import forcefield_rdk
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from simtk import unit
import parmed
from simtk import openmm
from pdbfixer import PDBFixer # for solvating
from openforcefield.typing.engines.smirnoff import forcefield_rdk
import sys

def create_system_from_molecule_rdk(forcefield, mol, verbose=False):
	"""
	Generate a System from the given OEMol and SMIRNOFF forcefield, return the resulting System.
	Parameters
	----------
	forcefield : ForceField
	    SMIRNOFF forcefield
	mol : RDKit molecule
	    Molecule to test (must have coordinates)
	Returns
	----------
	topology : OpenMM Topology
	system : OpenMM System
	positions : initial atomic positions (OpenMM)
	"""
	# Create system
	topology = utils.generateTopologyFromRDKMol(mol)
	system = forcefield.createSystem(topology, [mol], verbose=verbose)
	# Get positions
	coordinates = mol.GetConformer().GetPositions()
	natoms = len(coordinates)
	positions = np.zeros([natoms,3], np.float32)
	for index in range(natoms):
	    (x,y,z) = coordinates[index]
	    positions[index,0] = x
	    positions[index,1] = y
	    positions[index,2] = z
	positions = unit.Quantity(positions, unit.angstroms)
	return topology, system, positions


def generateSMIRNOFFSystemRDK(molecule):
	"""
	Given an RDKit molecule, create an OpenMM System and use to
	generate a ParmEd structure using the SMIRNOFF forcefield parameters.
	"""
	from openforcefield.typing.engines.smirnoff import forcefield_rdk
	from openforcefield.typing.engines.smirnoff.forcefield_utils import create_system_from_molecule
	#ff = utils.get_data_filename('forcefield/smirnoff99Frosst.ffxml')
	ff = 'smirnoff99Frosst.offxml'
	with open(ff) as ffxml:
	    mol_ff = forcefield_rdk.ForceField(ffxml)
	#TODO : integrate charges
	charged_molecule = molecule
	mol_top, mol_sys, mol_pos = create_system_from_molecule_rdk(mol_ff, charged_molecule)
	return mol_top, mol_sys, mol_pos


def generateSMIRNOFFStructureRDK(molecule):
	"""
	Given an RDKit molecule, create an OpenMM System and use to
	generate a ParmEd structure using the SMIRNOFF forcefield parameters.
	"""
	from openforcefield.typing.engines.smirnoff import forcefield_rdk
	from openforcefield.typing.engines.smirnoff.forcefield_utils import create_system_from_molecule
	#ff = utils.get_data_filename('forcefield/smirnoff99Frosst.ffxml')
	ff = 'smirnoff99Frosst.offxml'
	with open(ff) as ffxml:
	    mol_ff = forcefield_rdk.ForceField(ffxml)
	#TODO : integrate charges
	charged_molecule = molecule
	mol_top, mol_sys, mol_pos = create_system_from_molecule_rdk(mol_ff, charged_molecule)
	molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
	return molecule_structure



def parametrize_system(ligand_file,protein_file):
"""
	build new openmm system from ligand mol2 file and protein pdb file (pdb v3 format)
"""	
	lig_rdk=Chem.MolFromMol2File(ligand_file,sanitize=False)


	#m=Chem.MolFromSmiles('CC\C(=C(/c1ccc(O)cc1)c1ccc(OCCN(C)C)cc1)c1ccccc1',sanitize=True)
	#lig_rdk=Chem.AddHs(m)
	#AllChem.EmbedMolecule(lig_rdk, AllChem.ETKDG())
	lig_rdk.SetProp('_Name','LIG')

	ligand_pmd=generateSMIRNOFFStructureRDK(lig_rdk)


	print("fixing protein")
	protein=parmed.load_file(protein_file)
	protein.write_pdb("fixed.pdb")


	print("loading system")
	protein=parmed.load_file("fixed.pdb")
	protein = protein["!(:HOH,NA,CL)"] #remove ions and water

	forcefield = app.ForceField('amber99sb.xml')

	protein_system = forcefield.createSystem(protein.topology)

	protein_pmd = parmed.openmm.load_topology(protein.topology, protein_system, protein.positions)

	#pdb = app.PDBFile(protein_file) #openmm read-in pdb
	#forcefield = app.ForceField('amber99sb.xml', 'tip3p.xml')   #you can of course use other forcefields.xml files from openmm
	#protein_system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*unit.nanometer, constraints=app.HBonds)
	#protein_pmd = parmed.openmm.load_topology(pdb.topology, protein_system, pdb.positions)

	protein_pmd.save("protein_prepared.pdb",overwrite=True)
	#protein_pmd.save("protein_prepared.prmtop", overwrite=True)



	prot_lig_pmd = protein_pmd + ligand_pmd

	prot_lig_pmd.save("complex.pdb",overwrite=True)
	#prot_lig_pmd.save("complex.prmtop", overwrite=True)


	print("Solvation")
	fixer = PDBFixer("complex.pdb")


	#0.1 in Vec3 because box_size is in angstrom and fixer uses nanometer
	#scaling factor to somehow ensure no interaction with periodic image 
	scaling_factor = 1.0
	box_size=50
	fixer.addSolvent(scaling_factor * box_size *openmm.Vec3(0.1, 0.1, 0.1), positiveIon='Na+', negativeIon='Cl-', ionicStrength=0.1*unit.molar)

	app.PDBFile.writeFile(fixer.topology, fixer.positions, open('complex_solvated.pdb', 'w'))

	print("Solvation done")
	print("Parametrizing ions")

	complex = parmed.load_file('complex_solvated.pdb')

	ions = complex["(:NA,CL)"]


	forcefield = app.ForceField('amber99sb.xml')
	ions_system = forcefield.createSystem(ions.topology)
	# print(solvent_system.getNumConstraints())

	ions_pmd = parmed.openmm.load_topology(ions.topology, ions_system, ions.positions)
	print("Parametrizing ions done")


	prot_lig_ion_pmd = protein_pmd + ligand_pmd + ions_pmd

	#prot_lig_pmd.save("complex_ions.pdb",overwrite=True)
	#prot_lig_pmd.save("complex_ions.prmtop", overwrite=True)

	print("Parametrizing solvent")

	solvent = complex["(:HOH)"]
	num_solvent = len(solvent.residues)


	solvent_pmd = parmed.load_file("../../water/water.prmtop")
	solvent_pmd *= num_solvent
	solvent_pmd.positions = solvent.positions
	print("Parametrizing solvent done")


	print("merge structures")

	combined_pmd = protein_pmd + ligand_pmd + ions_pmd + solvent_pmd
	#combined_pmd.write_pdb("system_complex.pdb",renumber=False)
	combined_pmd.save("system_complex.prmtop", overwrite=True)
	print("merge done")
	combined_pmd.box_vectors = complex.box_vectors
	

	print("writing pickle")
	import pickle
	pickle_out = open("complex_system.pickle","wb")
	pickle.dump([combined_pmd], pickle_out)
	pickle_out.close()
	
	print("system successfully parametrized")



if __name__ == "__main__":
    # execute only if run as a script

	if len(sys.argv)!=3 : 
		sys.exit("USAGE : python 01_parametrize.py protein.pdb ligand.mol2\n")
	ligand_file=sys.argv[2]
	protein_file=sys.argv[1]

	parametrize_system(ligand_file, protein_file)



