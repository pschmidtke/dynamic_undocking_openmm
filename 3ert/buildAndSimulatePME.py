from simtk.openmm import app
import openforcefield.utils as utils
from openforcefield.typing.engines.smirnoff import forcefield_rdk
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from simtk import unit
import parmed
from simtk import openmm
from simtk.openmm.app import PDBFile
from pdbfixer import PDBFixer 
#import packmol

from openforcefield.typing.engines.smirnoff import forcefield_rdk

platform = openmm.Platform.getPlatformByName('OpenCL')


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

def generateSolventStructure(solventpdb, solvent_forcefield='tip3p.xml'):
	"""
	Given an OpenMM PDBFile, create the OpenMM System of the protein and
	then generate the parametrized ParmEd Structure of the protein.
	Parameters
	----------
	proteinpdb : openmm.app.PDBFile object,
	    Loaded PDBFile object of the protein.
	protein_forcefield : xml file, default='amber99sbildn.xml'
	    Forcefield parameters for protein
	solvent_forcefield : xml file, default='tip3p.xml'
	    Forcefield parameters for solvent
	Returns
	-------
	solv_structure : parmed.structure.Structure
	    The parameterized Structure of the protein with solvent molecules. (No ligand).
	"""
	#Generate protein Structure object
	from openforcefield.typing.engines.smirnoff import forcefield_rdk
	from openforcefield.typing.engines.smirnoff.forcefield_utils import create_system_from_molecule
	ff = 'tip3p.offxml'
	with open(ff) as ffxml:
		mol_ff = forcefield_rdk.ForceField(ffxml)
	solvent_system = mol_ff.createSystem( solventpdb.topology ,[])
	solvent_structure = parmed.openmm.load_topology(solventpdb.topology,solvent_system,xyz=solventpdb.positions)
	return solvent_structure




lig_rdk=Chem.MolFromMol2File('ligand.mol2',sanitize=True)
lig_rdk.SetProp('_Name','lig')
print("parametrizing ligand")
molecule_structure=generateSMIRNOFFStructureRDK(lig_rdk)

molecule_structure.save("ligand.pdb",overwrite=True)
molecule_structure.save("ligand.prmtop", overwrite=True)

#print("adding solvent to structure")
#fixer=PDBFixer("3ert.pdb")
#fixer.addSolvent(padding=unit.Quantity( 1.0, unit.angstroms), ionicStrength=unit.Quantity( 20, unit.millimolar))

#PDBFile.writeFile(fixer.topology, fixer.positions, open('complex_solvated.pdb', 'w'))
tmp=Chem.MolFromSmiles('O',sanitize=True)
water_rdk=Chem.AddHs(tmp)
water_rdk.SetProp('_Name','wat')
print('parametrizing solvent')
solventpdb = app.PDBFile('solvent.pdb')
solventmin=solventpdb.getPositions().min()
solventmax=solventpdb.getPositions().max()
from openforcefield.typing.engines import smirnoff
ff = smirnoff.forcefield_rdk.ForceField('smirnoff99Frosst.offxml', 'tip3p.offxml')
system_solvated = ff.createSystem(solventpdb.topology, molecules=[water_rdk],nonbondedMethod=smirnoff.PME, nonbondedCutoff=1.1*unit.nanometer,ewaldErrorTolerance=1e-4) 



solvent_structure = generateSolventStructure(solventpdb)
print(solvent_structure)

print("building protein system")
proteinpdb = app.PDBFile('3ert.pdb')
#parmedpdb=parmed.load_file('1uyd.pdb')
protein_structure = utils.generateProteinStructure(proteinpdb)
structure=protein_structure
#print("merging with ligand")
#structure = utils.mergeStructure(protein_structure, molecule_structure)

structure.save("system.pdb",overwrite=True)
#structure.save("system.prmtop", overwrite=True)



print("creating simulation system")
#system=protein_structure.createSystem(constraints=None, rigidWater=None,flexibleConstraints=None)
#system=structure.createSystem( nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)
system=structure.createSystem(flexibleConstraint=False)

integrator = openmm.LangevinIntegrator(300*unit.kelvin,91/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(structure.topology, system, integrator,platform)
simulation.context.setPositions(structure.positions)
print("starting minimization")
simulation.minimizeEnergy()
simulation.reporters.append(app.DCDReporter("output.dcd",1000))
print("starting simulation")
simulation.step(100000000)








