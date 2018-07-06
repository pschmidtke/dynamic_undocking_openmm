import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
#!#from mdtraj.reporters import HDF5Reporter 
import numpy as np
import sys
import pickle
import duck




if len(sys.argv)!=7 :
  sys.exit("Usage 04_smd.py temperature in.chk out.csv out.dat out.pdb out.dcd")




temperature = float(sys.argv[1])*u.kelvin
checkpoint_in_file = sys.argv[2]
csv_out_file = sys.argv[3]
dat_out_file = sys.argv[4]
pdb_out_file = sys.argv[5]
traj_out_file = sys.argv[6]

#keyInteraction_ind_mol2 = [453, 1329]
#keyInteraction = [453, 1329] #[keyInteraction_ind_mol2[0]-1, keyInteraction_ind_mol2[1]-1]

spring_k = 50 * u.kilocalorie/(u.mole * u.angstrom * u.angstrom)
dist_in = 2.5 * u.angstrom # in angstrom 
dist_fin = 5.0 * u.angstrom # in angstrom

steps_per_move = 200

velocity = 0.00001 * u.angstrom

# Platform definition

platform = mm.Platform_getPlatformByName("OpenCL")
platformProperties = {}
platformProperties['OpenCLPrecision'] = 'mixed'
#platformProperties['CudaDeviceIndex'] = '0'

print("loading pickle")
pickle_in=open('complex_system.pickle', 'rb')
combined_pmd = pickle.load(pickle_in)[0]
print(combined_pmd)
pickle_in.close()



resnumber="29"
atomname="O"
distance="3.0"

keyInteraction=duck.getAtomSerialFromAmberMask(combined_pmd,resnumber,atomname,distance)
print(keyInteraction)


# Get indexes of heavy atoms in chunk

Chunk_Heavy_Atoms = duck.getCAAtomsInSystem(combined_pmd)

# Setting System

system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)


# Apply force on all havy atoms of chunk

duck.applyHarmonicPositionalRestraints(system, 0.5, combined_pmd.positions, Chunk_Heavy_Atoms)


# Integrator

integrator = mm.LangevinIntegrator(temperature, 4/u.picosecond, 0.002*u.picosecond)


# Setting Simulation object and loading the checkpoint

simulation = app.Simulation(combined_pmd.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Get positions

positions = simulation.context.getState(getPositions=True).getPositions()


# SMD force definition

pullforce = mm.CustomExternalForce('k_sp*0.5*(R-R0)^2; \
                                   R = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);')

pullforce.addPerParticleParameter('k_sp')

pullforce.addGlobalParameter("x0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("y0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("z0", 0.0 * u.nanometer)
pullforce.addGlobalParameter("R0", 0.0 * u.nanometer)

pullforce.addParticle(keyInteraction[1], [spring_k])

system.addForce(pullforce)


# Redefine integrator and simulation, and load checkpoint with new-updated system

integrator = mm.LangevinIntegrator(temperature, 4/u.picosecond, 0.002*u.picosecond)
simulation = app.Simulation(combined_pmd.topology, system, integrator, platform, platformProperties)
simulation.loadCheckpoint(checkpoint_in_file)


# Initializing energy

work_val_old = u.Quantity(value=0, unit=u.kilocalorie/u.mole)


# Number of big steps and pull distance 

steps = int(round((dist_fin  - dist_in) / velocity) / steps_per_move)
pull_distance = velocity * steps_per_move


# Reporters and duck.dat file

simulation.reporters.append(app.StateDataReporter(csv_out_file, steps_per_move, step=True, time=True, totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, progress=True, totalSteps = steps_per_move * steps, speed=True))
simulation.reporters.append(app.DCDReporter(traj_out_file, 1000))
f=open(dat_out_file,'w')


# Production in N steps with the update every 200 steps (2 pm) 

for i in range(steps):

    # Get current state tu update the system

    state = simulation.context.getState(getPositions=True)
    pos_keyInt = state.getPositions()
    keyInteraction_pos = [pos_keyInt[keyInteraction[0]], pos_keyInt[keyInteraction[1]]]

    # Get radius of starting point and end point

    R_val = (dist_in + float(i+1) * pull_distance)
    R_val_start = (dist_in + float(i) * pull_distance)

    # Get distance of main interaction

    keyInteraction_dist = np.linalg.norm(keyInteraction_pos[0]-keyInteraction_pos[1])
    print(keyInteraction_dist)
    # Updated system

    simulation.context.setParameter('x0', keyInteraction_pos[0][0])
    simulation.context.setParameter('y0', keyInteraction_pos[0][1])
    simulation.context.setParameter('z0', keyInteraction_pos[0][2])
    simulation.context.setParameter('R0', R_val)

    # Calculate force F = -k * x

    force_val = -spring_k * (keyInteraction_dist - R_val)

    # Make step

    simulation.step(steps_per_move)

    # Calculate work for difference in potential energy in tranzition
    # W = EK_end - EK_start 
    # EK = 0.5 * k * x^2

    spr_energy_end = 0.5 * spring_k * (keyInteraction_dist - R_val)**2
    spr_energy_start = 0.5 * spring_k * (keyInteraction_dist - R_val_start)**2
    
    work_val = work_val_old + spr_energy_end - spr_energy_start

    work_val_old = work_val

    # Write duck.dat file 

    f.write(str(i)+' '+str(R_val)+' '+str(keyInteraction_dist)+' '+str(force_val)+' '+str(work_val)+'\n')
    
    
f.close()

# Save state in PDB file

positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, 'w'))
