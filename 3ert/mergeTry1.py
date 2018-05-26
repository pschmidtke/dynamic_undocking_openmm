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


lig_rdk=Chem.MolFromMol2File('ligand.mol2',sanitize=True)


#m=Chem.MolFromSmiles('CC\C(=C(/c1ccc(O)cc1)c1ccc(OCCN(C)C)cc1)c1ccccc1',sanitize=True)
#lig_rdk=Chem.AddHs(m)
#AllChem.EmbedMolecule(lig_rdk, AllChem.ETKDG())
lig_rdk.SetProp('_Name','mymol')

molecule_structure=generateSMIRNOFFStructureRDK(lig_rdk)




molecule_structure.save("ligand.pdb",overwrite=True)
molecule_structure.save("ligand.prmtop", overwrite=True)



#fixer = PDBFixer( "3ert.pdb")


#fixer.addSolvent(padding=unit.Quantity( 10.0, unit.angstroms), ionicStrength=unit.Quantity( 20, unit.millimolar))
#print(dir(fixer))
#f=open("system_solvated.pdb","w")
#app.PDBFile.writeFile(fixer.topology, fixer.positions, f, True)
#f.close()

proteinpdb = app.PDBFile('system_solvated.pdb')

#parmedpdb=parmed.load_file('1uyd.pdb')

protein_structure = utils.generateProteinStructure(proteinpdb)

structure = utils.mergeStructure(protein_structure, molecule_structure)
structure.save("system.pdb",overwrite=True)
structure.save("system.prmtop", overwrite=True)



#system=protein_structure.createSystem(constraints=None, rigidWater=None,flexibleConstraints=None)
system=structure.createSystem(constraints=None)

integrator = openmm.LangevinIntegrator(300*unit.kelvin,1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(structure.topology, system, integrator)
simulation.context.setPositions(structure.positions)	
simulation.minimizeEnergy()
simulation.reporters.append(app.DCDReporter("output.dcd",10))
simulation.step(10000)


mol_top, mol_sys, mol_pos=generateSMIRNOFFSystemRDK(lig_rdk)
integrator = openmm.LangevinIntegrator(300*unit.kelvin,1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(mol_top, mol_sys, integrator)
simulation.context.setPositions(mol_pos)
simulation.minimizeEnergy()
simulation.reporters.append(app.DCDReporter("output2.dcd",10))
simulation.step(10000)








