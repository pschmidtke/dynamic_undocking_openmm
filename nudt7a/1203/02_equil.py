import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import pickle
import duck #set up pythonpath first to include duck.py to get all helper functions

###############################################################################

###############################################################################
#################### Sym parameters ###########################################
###############################################################################

keyInteraction_ind_mol2 = [453, 1329]
keyInteraction = [keyInteraction_ind_mol2[0]-1, keyInteraction_ind_mol2[1]-1]

# Platform definition

platform = mm.Platform_getPlatformByName("CPU")
platformProperties = {}
#platformProperties['OpenCLPrecision'] = 'mixed'
#platformProperties['CudaDeviceIndex'] = '0'


# Read files

print("loading pickle")
pickle_in=open('complex_system.pickle', 'rb')
combined_pmd = pickle.load(pickle_in)[0]
print(dir(combined_pmd))
pickle_in.close()

##################
##################
#  Minimisation  #
##################
##################

print('Minimising...')

# Define system

system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom)

# Get indexes of heavy atoms in chunk

Chunk_Heavy_Atoms = duck.getHeavyAtomsInSystem(combined_pmd)

# Apply force on all havy atoms of chunk

duck.applyHarmonicPositionalRestraints(system, 1.0, combined_pmd.positions, Chunk_Heavy_Atoms)

# Integrator

integrator = mm.VerletIntegrator(1*u.femtosecond)

# Define Simulation

simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
simulation.context.setPositions(combined_pmd.positions)

# Minimizing energy

simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions

positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))


##########################
##########################
# Equlibration - heating #
##########################
##########################
#new minimised positions, however using old restraints

# Define new system

system = combined_pmd.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=9*u.angstrom, constraints=app.HBonds, hydrogenMass=None)

# Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance

duck.applyHarmonicPositionalRestraints(system, 1.0, combined_pmd.positions, Chunk_Heavy_Atoms)
duck.applyLigandChunkRestraint(system, 1.0, 10.0, 2*u.angstrom, 3*u.angstrom, 4*u.angstrom, keyInteraction)

# Intergator 

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

# Define Simulation

simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
simulation.context.setPositions(positions) #changing coordintes to minimized   

# Reporters

simulation.reporters.append(app.StateDataReporter("heating.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Heating the system
print("Heating ... ")

simulation.step(50000) # 0.1 ns 

# Save the positions and velocities

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

app.PDBFile.writeFile(simulation.topology, positions, open('heating_final.pdb', 'w'))

#clear reporters

simulation.reporters = []

##########################
##########################
# Equlibration - density #
##########################
##########################

# Add barostat to the system

system.addForce(mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin))

# Integrator 

integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

# Define simulation

simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

# Reporters

simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Correcting the density
print("Correcting density")

simulation.step(50000) # 0.1 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))

#saving simulation stage
simulation.saveCheckpoint('equil.chk')
