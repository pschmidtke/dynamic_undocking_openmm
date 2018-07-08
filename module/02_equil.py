import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import pickle
import duck #set up pythonpath first to include duck.py to get all helper functions
import sys

###############################################################################

###############################################################################
#################### Sym parameters ###########################################
###############################################################################

if len(sys.argv)!=5 :
  sys.exit("Usage 02_equil.py resnumber atomname distance hbond(1/0)")


resnumber=sys.argv[1]
atomname=sys.argv[2]
distance=sys.argv[3]
hbond=sys.argv[4]



# Platform definition

platform = mm.Platform_getPlatformByName("CPU")


platformProperties = {}
#platformProperties['OpenCLPrecision'] = 'mixed'
#platformProperties['UseCpuPme'] = 'false'
#platformProperties['DisablePmeStream'] = 'true'
#platformProperties['CudaDeviceIndex'] = '0'


# Read files

print("loading pickle")
pickle_in=open('complex_system.pickle', 'rb')
combined_pmd = pickle.load(pickle_in)[0]
combined_pmd.symmetry=None
pickle_in.close()


if hbond=="1":
	keyInteraction=duck.getAtomSerialFromAmberMaskHbond(combined_pmd,resnumber,atomname,distance)
else :
	keyInteraction=duck.getAtomSerialFromAmberMask(combined_pmd,resnumber,atomname,distance)
	
print(keyInteraction)
#print("selecting ligand atom index ")
#lig_sel=parmed.amber.mask.AmberMask(combined_pmd,"(((:29)&(@O))<@3.0) & (:LIG & (@N=|@O=))")
##lig_idx=[x for x in lig_sel.Selected() ]
#print("selecting protein atom index ")
#rec_sel=parmed.amber.mask.AmberMask(combined_pmd,"((:29)&(@O))")
#rec_idx=[x for x in rec_sel.Selected() ]

#if(len(lig_idx)>1 or len(rec_idx)>1) :
#	sys.exit("The reaction coordinate selection failed here, please consider setting it by hand")

#keyInteraction = [lig_idx[0], rec_idx[0]]
#sys.exit()

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

print(simulation.context.getPlatform().getName())
for key in simulation.context.getPlatform().getPropertyNames():
    print(key, simulation.context.getPlatform().getPropertyValue(simulation.context, key))



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
#system = combined_pmd.createSystem(nonbondedMethod=app.PME)




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
simulation.reporters.append(app.DCDReporter("heating.dcd", 1000))

# Heating the system
print("Heating ... ")

simulation.step(5000) # 0.01 ns 

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

simulation = duck.setUpNPTEquilibration(system,combined_pmd,platform,platformProperties,positions,velocities)
# Reporters

simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Correcting the density
print("Correcting density CPU")

simulation.step(10000) # 0.01 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))






platform = mm.Platform_getPlatformByName("OpenCL")


platformProperties = {}
platformProperties['OpenCLPrecision'] = 'mixed'
#platformProperties['UseCpuPme'] = 'false'
#platformProperties['DisablePmeStream'] = 'true'
#platformProperties['CudaDeviceIndex'] = '0'



simulation.reporters = []

##########################
##########################
# Equlibration - density #
##########################
##########################

# Add barostat to the system

#system.addForce(mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin))

# Integrator 

#integrator = mm.LangevinIntegrator(300*u.kelvin, 4/u.picosecond, 0.002*u.picosecond)

# Define simulation

#simulation = app.Simulation(combined_pmd.topology, system, integrator, platform,platformProperties)
#simulation.context.setPositions(positions)
#simulation.context.setVelocities(velocities)


simulation = duck.setUpNPTEquilibration(system,combined_pmd,platform,platformProperties,positions,velocities)

# Reporters

simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=50000))

# Correcting the density
print("Correcting density OpenCL")

simulation.step(5000) # 0.01 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
app.PDBFile.writeFile(simulation.topology, positions, open('density_final.pdb', 'w'))


#saving simulation stage
simulation.saveCheckpoint('equil.chk')
