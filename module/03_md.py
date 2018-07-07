import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
#!#from mdtraj.reporters import HDF5Reporter
import sys
import pickle
import duck

if len(sys.argv)!=9 :
  sys.exit("Usage 03_md.py in.chk out.chk out.csv out.pdb resnumber atomname distance hbond(1/0)")

print("loading pickle")
pickle_in=open('complex_system.pickle', 'rb')
combined_pmd = pickle.load(pickle_in)[0]
print(dir(combined_pmd))
pickle_in.close()



#resnumber="29"
#atomname="O"
#distance="3.0"


checkpoint_in_file = sys.argv[1]
checkpoint_out_file = sys.argv[2]
csv_out_file = sys.argv[3]
pdb_out_file = sys.argv[4]
resnumber=sys.argv[5]
atomname=sys.argv[6]
distance=sys.argv[7]
hbond=sys.argv[8]

traj_out_file = "traj.out" #sys.argv[5]
if hbond=="1":
	keyInteraction=duck.getAtomSerialFromAmberMaskHbond(combined_pmd,resnumber,atomname,distance)
else :
	keyInteraction=duck.getAtomSerialFromAmberMask(combined_pmd,resnumber,atomname,distance)
	
print(keyInteraction)

MD_len = 0.05 * u.nanosecond
sim_steps = round(MD_len / (0.002 * u.picosecond))



# Platform definition

platform = mm.Platform_getPlatformByName("OpenCL")
platformProperties = {}
platformProperties['OpenCLPrecision'] = 'mixed'




# Get indexes of heavy atoms in chunk

Chunk_Heavy_Atoms = duck.getHeavyAtomsInSystem(combined_pmd)

        
# Setting System

system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)


# Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance

duck.applyHarmonicPositionalRestraints(system, 0.1, combined_pmd.positions, Chunk_Heavy_Atoms)
duck.applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)


# Integrator

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)


# Setting Simulation object and loading the checkpoint

simulation = app.Simulation(combined_pmd.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Simulation reporters

simulation.reporters.append(app.StateDataReporter(csv_out_file, 2000, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps=250000, speed=True))
#!#simulation.reporters.append(HDF5Reporter(traj_out_file, XXX))
simulation.reporters.append(app.DCDReporter("md.dcd", 1000))


# Production

simulation.step(sim_steps)


# Save state in checkpoint file and save coordinates in PDB file
            
state = simulation.context.getState(getPositions=True, getVelocities=True)
positions = state.getPositions()
velocities = state.getVelocities()

app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, 'w'))

simulation.saveCheckpoint(checkpoint_out_file)
