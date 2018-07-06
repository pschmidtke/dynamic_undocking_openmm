import simtk.openmm as mm
import simtk.unit as u
import parmed

#!#from mdtraj.reporters import HDF5Reporter 



def applyHarmonicPositionalRestraints(system, forceConstantInKcalPerMolePerAngSquared,
                                      positions, indexOfAtomsToBeModified):
    """ This is essentially mimicking AMBER's restraint_wt"""

    forceConstant = u.Quantity(value=forceConstantInKcalPerMolePerAngSquared,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))

    force = mm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")


    force.addGlobalParameter("k",
       forceConstant.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))

    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for i in indexOfAtomsToBeModified:
        force.addParticle(i, positions[i])

    system.addForce(force)


def applyLigandChunkRestraint(system, k2, k3, R2, R3, R4, indexOfAffectedAtoms):
    
    forceConstant_k2 = u.Quantity(value=k2,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))
    forceConstant_k3 = u.Quantity(value=k3,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))
    
    restraint_force = mm.CustomBondForce('step(R2 - r) * f1 + step(r - R3) * select( step(r - R4), f3, f2);'
                                         'f1 = k2 * (r - R2)^2;'
                                         'f2 = k3 * (r - R3)^2;'
                                         'f3 = k3 * (R4 - R3) * (2 * r - R4 - R3)')
    
    restraint_force.addGlobalParameter("k2",
       forceConstant_k2.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))
    restraint_force.addGlobalParameter("k3",
       forceConstant_k3.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))
    
    restraint_force.addGlobalParameter("R2", R2) 
    restraint_force.addGlobalParameter("R3", R3)
    restraint_force.addGlobalParameter("R4", R4)
    
    restraint_force.addBond(indexOfAffectedAtoms[0], indexOfAffectedAtoms[1])
    
    system.addForce(restraint_force)


def loadPickle(fname='complex_system.pickle'):
    """ load parametrized system from pickle to continue to set up simulations
        fname = optional parameter which file to read
        careful a pickle contains also the openmm platform parameters. 
        So pickles created with CPU can only be used for CPU based simulation. Same for CUDA and OpenCL.
    """
    pickle_in=open(fname, 'rb')
    combined_pmd = pickle.load(pickle_in)[0]
    pickle_in.close()
    return(combined_pmd)

def getHeavyAtomsInSystem(combined_pmd):
    Chunk_Heavy_Atoms=[]
    for atom_i in combined_pmd.topology.atoms():
        if atom_i.residue.name not in ('HOH', 'WAT', 'IP3', 'LIG', 'Cs+', 'K+', 'Rb+', 'Li+', 'Na+', 'IP', 'Cl-', 'IM', 'IB') and atom_i.name[0] != 'H':
            Chunk_Heavy_Atoms.append(atom_i.index)
    return(Chunk_Heavy_Atoms)

def getCAAtomsInSystem(combined_pmd):
    Chunk_Heavy_Atoms=[]
    for atom_i in combined_pmd.topology.atoms():
        if atom_i.residue.name not in ('HOH', 'WAT', 'IP3', 'LIG', 'Cs+', 'K+', 'Rb+', 'Li+', 'Na+', 'IP', 'Cl-', 'IM', 'IB') and atom_i.name[0] != 'CA':
            Chunk_Heavy_Atoms.append(atom_i.index)
    return(Chunk_Heavy_Atoms)


def getAtomSerialFromAmberMask(combined_pmd,resnumber,atomname,distance,ligresname="LIG"):
    lig_sel=parmed.amber.mask.AmberMask(combined_pmd,"(((:"+resnumber+")&(@"+atomname+"))<@3.0) & (:"+ligresname+" & (@N=|@O=))")
    lig_idx=[x for x in lig_sel.Selected() ]
    rec_sel=parmed.amber.mask.AmberMask(combined_pmd,"((:"+resnumber+")&(@"+atomname+"))")
    rec_idx=[x for x in rec_sel.Selected() ]


    if(len(lig_idx)>1 or len(rec_idx)>1) :
      sys.exit("The reaction coordinate selection failed here, please consider setting it by hand")
    return(rec_idx[0],lig_idx[0])

