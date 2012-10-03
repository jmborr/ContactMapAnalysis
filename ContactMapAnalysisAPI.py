import MDAnalysis as MDA
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.KDTree.NeighborSearch import AtomNeighborSearch
from scipy.sparse import csr_matrix
import ContactMapAnalysis as CMA
import numpy
import sys
from pdb import set_trace as trace

"""
Desired properties:
(1) parser defining how to evaluate a contact
(2) iterate over the trajectory frames
(3) compute the contact map with the given parser and frame, returning a crs_ContactMap object
"""

def AtomsListPerResidue(atoms,shift=True):
    """Returns a list whose elements are also lists. Each element contains the indices                                  
    of the atoms contained in 'atoms' that belong to a particular residue.                                              
    Note: not all atoms belonging to a particular residue may be contained in 'atoms'                                   
    shift=True, by default, numbering of atoms and residues is lost.                                                    
          Example: if 'atoms' contains atoms with numbers 3,5,8,10,12, the numbers become 0,1,2,3,4.                    
                   Similary, if 'atoms' contains atoms in residues 5,9,12, the number become 0,1,2                      
    """
    ResDict={} #as many entries as residues in atoms                                                               
    for resNum in atoms.resnums():
        ResDict[resNum]=[]
    iat=0 #we loose the atom number convention in the PSF file                                                          
    for atom in atoms:
        ResDict[ atom.resid ].append(iat)
        iat+=1
    ResList=[]
    for key in sorted(ResDict.iterkeys()): ResList.append(ResDict[key])
    return ResList    
    
class ContactMapProtocol:
    """protocol to generate the contact map"""
    def __init__(self,**kwargs):
        """Attributes:
          cutOff, cut off defining the contact, in Angstroms
          rowSelection, atom selection to be the row index of the Contact map
          colSelection, atom selection to be the column index of the Contact map
          rowGroup, atomGroup object, from applying selectAtoms to rowSelection
          colGroup, atomGroup object, from applying selectAtoms to colSelection
          self.shift, set to 1 if the first atomic index in Universe is zero
          byres=False, if set to True, contact map between residues instead of between atoms
          reducer=numpy.any, algorithm to define the contact between two residues based 
                  on the atomic contact map between the two residues
                  Example: if residue i made up of atoms 0,1,2 and residue j made up of atoms 4,5
                           then C[[0,1,2]][:,[4,5]] is the desired C_ij submatrix.
                  numpy.any( C_ij ) is evaluated by default, but other reducers can be supplied
          filter: python function to filter non-wanted contacts (e.g self-contacts). Takes a csr_matrix as argument
          isInitialized=False, initialized() hast not been called
        """
        self.isInitialized=False
        self.byres=False
        self.reducer=numpy.any
        """
        #One example of filtering contacts
        def filter(csr):
            # filter out contacts that are too close along the sequence
            # csr: contact map, a csr_matrix object
            cutOff=1
            i=0
            for irow in range( len(csr.indptr) - 1 ):
                for icol in csr.indices[ csr.indptr[irow]:csr.indptr[irow+1]]:
                    if abs(irow-icol)<=cutOff: csr.data[i] = 0
                    i += 1
            csr.eliminate_zeros()
            return csr
        """
        self.filter = None
        self.__dict__.update(kwargs)
        
    def initialize(self,Universe):
        if self.isInitialized: return None
        try:
            self.rowGroup = Universe.selectAtoms(self.rowSelection)
            self.colGroup = Universe.selectAtoms(self.colSelection)
            self.shift=0
            if Universe.atoms[0].number==0: self.shift=1
        except AttributeError:
            print 'rowSelection and/or colSelection undefined in protocol\n'
        if self.byres:
            """generate atom number to residue number conversion array"""
            self.numberOfResidues=Universe.residues.numberOfResidues()
            tmp={}
            for atom in Universe.atoms: tmp[atom.number]=atom.resid
            N=max(tmp.keys())+self.shift
            self.atom2res=numpy.zeros(N)
            for (atnum,resid) in tmp.items(): self.atom2res[atnum]=resid
            
def GenerateContactMap(TimeStep,protocol):
    """Read a frame and compute the contactMap between two groups of atom given a ContactMapProtocol
    NOTE: if protocol.byres=True, then we report the number of atomic contacts for a given
          residue pair in contact
    Numbering of atoms and residues is the one given by the PSF file"""
    rowNS = AtomNeighborSearch(protocol.rowGroup) #KDE-tree object for rowGroup
    colNS = AtomNeighborSearch(protocol.colGroup) #KDE-tree object for colGroup
    rowClose = rowNS.search_list(protocol.colGroup,protocol.cutOff) #atoms of rowGroup in contact with colGroup
    if not rowClose: #no contacts
        shape=(TimeStep.numatoms,TimeStep.numatoms) #all atoms in the system
        if protocol.byres: shape=(protocol.numberOfResidues,protocol.numberOfResidues)
        csr=csr_matrix( ([],([],[])), shape, dtype='int32')  #empty map
        cma=CMA.csr_ContactMap(csr)
        cma.setDistanceCutOff(protocol.cutOff)
        return cma
    colClose = colNS.search_list(rowClose,protocol.cutOff) #atoms of colGroup in contact with rowClose/rowGroup
    dd=distance_array(rowClose.coordinates(),colClose.coordinates(),TimeStep.dimensions[0:3])
    (rowIndices,colIndices) = numpy.where(dd<protocol.cutOff)
    rowIndices=rowClose.indices()[rowIndices]
    colIndices=colClose.indices()[colIndices]  
    if protocol.byres:
        #switch from atomic indices to residue numbers
        rowIndices=protocol.atom2res[rowIndices]
        colIndices=protocol.atom2res[colIndices]
        shape=(protocol.numberOfResidues,protocol.numberOfResidues)
    else:
        #Take into account if the first atomic index was zero
        rowIndices += protocol.shift
        colIndices += protocol.shift
        shape=(TimeStep.numatoms,TimeStep.numatoms) #all atoms in the system
    #create sparse matrix for the contact map. 
    data=numpy.ones( len(rowIndices) ) #just signaling a contact map
    csr=csr_matrix( (data, (rowIndices,colIndices)), shape, dtype='int32')
    if protocol.byres: csr.data[:]=1 #overwrite number of atomic contacts per residue pair with 1 only
    #trace()
    if protocol.filter: csr=protocol.filter(csr) #Filtering
    cma=CMA.csr_ContactMap(csr)
    cma.setDistanceCutOff(protocol.cutOff)
    return cma
    
def GenerateContactMapList(PSFfile,trajfile,protocol):
    """Read a trajectory and compute the contactMap between two groups of atoms
    given a ContactMapProtocol"""
    u=MDA.Universe(PSFfile,trajfile)
    protocol.initialize(u)
    cml=CMA.ContactMapList()
    print 'Reading frame...', ; sys.stdout.flush()
    for ts in u.trajectory: # iterate through all frames
        cm=GenerateContactMap(ts,protocol)
        cml.append(cm)
        if not ts.frame%1000: print ts.frame,'..', ; sys.stdout.flush()
    print 'finished!'
    return cml
    
