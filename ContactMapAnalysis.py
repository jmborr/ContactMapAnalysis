import numpy
import h5py
import sys
from scipy.sparse import csr_matrix
from pdb import set_trace as trace
from copy import copy

def loadResidenceTimesListFromHDF5(Hdf5Group):
    rtl=ResidenceTimesList( Hdf5Group.attrs['n'])
    rtl.setDt( Hdf5Group.attrs['dt'] )
    for dataset in Hdf5Group:
        idx=int(dataset)
        rtl.data[idx]=list(Hdf5Group[dataset].value)
    return rtl
    
def loadResidenceTimesListFromFile(filename,fmt='HDF5'):
    """load list of contact maps from a file"""
    rtl = None
    extDict={'hd5':'HDF5','txt':'DAT','dat':'DAT'}
    if not fmt:
        ext=filename.split('.')[-1]
        fmt=extDict[ext]
    if fmt=='HDF5':
        f = h5py.File(filename,'r')
        rtl = loadResidenceTimesListFromHDF5(f['/'])
    else:
        print 'file format not understood or yet implemented'
    return rtl

class ResidenceTimesList:
    """a 2D array of residence time lists"""
    def __init__(self,n):
        self.n=n
        self.data=numpy.empty((n,),dtype=object)
        self.dt=1.0 #unit of time, in picoseconds
        for i in range(n): self.data[i]=[]
        pass
    def setDt(self,dt): self.dt = float(dt)  #set the time unit
    def applyFunc(self,myfunc):
        """apply myfunc to each element of ContactMapList"""
        vfunc=numpy.vectorize(myfunc)
        return vfunc(self)
    def append(self,cm,TimeCutOff):
        """append *positive* residence times from a 2D matrix 'cm'
        TimeCutOff: append times only bigger than TimeCutOff"""
        #silly implementation
        for i in range(self.n):
            A=cm[i].data
            self.data[i] += list(A[A>TimeCutOff])
    def saveToHDF5(self,Hdf5Group):
        """save the list to a particular group in the HDF5 file"""
        Hdf5Group.attrs['n']=self.n
        Hdf5Group.attrs['dt']=self.dt
        for idx in range(self.n):
            if self.data[idx]: #we cannot save empty lists
                Hdf5Group.create_dataset('%05d'%idx,data=self.data[idx])
    def saveToFile(self,filename,mode='w',fmt='HDF5'):
        """save to file in the selected format"""
        if fmt=='HDF5':
            f = h5py.File(filename,mode)
            self.saveToHDF5(f['/'])
        else:
            sys.stderr.write('ERROR: format '+fmt+' not yet implemented\n')

                    
class csr_ContactMap(csr_matrix):
    """a typedef from csr_matrix.
    The 'csr_' prefix is *required* for correct subclassing"""
    def __new__(cls,*kargs, **kwargs):
        obj=csr_matrix.__new__(cls, *kargs, **kwargs)
        return obj
    def setDistanceCutOff(self,co): self.cutOff=co
    def saveToHDF5(self,Hdf5Group):
        """save the contact map to a particular group in the HDF5 file"""
        Hdf5Group.attrs['cutOff']=self.cutOff
        Hdf5Group.attrs['shape']=self.shape
        A=self.tocoo() #cast to coo_matrix to retrieve row and col attributes
        #unfortunately, HDF5 doesn't allow saving empty map
        if A.data.any(): #map is not empty
            Hdf5Group.create_dataset('row',data=A.row)
            Hdf5Group.create_dataset('col',data=A.col)
            Hdf5Group.create_dataset('data',data=A.data)
    def saveToFile(self,filename,mode='a',fmt='HDF5',index=0):
        """save to file in the selected format"""
        if fmt=='HDF5':
            f = h5py.File(filename,mode)
            Hdf5Group=f.create_group('%05d'%index)
            Hdf5Group.attrs['index']=index
            self.saveToHDF5(Hdf5Group)
        else:
            sys.stderr.write('ERROR: format '+fmt+' not yet implemented\n')

def load_csr_ContactMapFromHDF5(Hdf5Group):
    row=[] ; col=[]; data=[]
    if 'data' in Hdf5Group.keys():
        row=Hdf5Group['row'].value
        col=Hdf5Group['col'].value
        data=Hdf5Group['data'].value      
    A=csr_matrix( (data, (row,col)), shape=Hdf5Group.attrs['shape'] )
    B=csr_ContactMap(A)
    B.setDistanceCutOff(Hdf5Group.attrs['cutOff'])
    return B

def load_csr_ContactMapFromFile(filename,fmt='HDF5',index=0):
    if fmt=='HDF5':
        f = h5py.File(filename,'r')
        Hdf5Group=f['%05d'%index]
        return load_csr_ContactMapFromHDF5(Hdf5Group)
    sys.stderr.write('ERROR: format '+fmt+' not yet implemented\n')
    return None
    
class ContactMapList(list):
    "a python list of csr_ContactMap objects"     
    def __new__(cls,*kwargs):
        obj=list.__new__(cls,*kwargs)
        try: 
            obj.shape=obj[0].shape
            obj.cutOff=obj[0].cutOff
        except:
            obj.shape=None
            obj.cutOff=None
        return obj
    def setShape(self,shape): self.shape=shape
    def setCutOff(self,cutOff): self.cutOff=cutOff
    def Occupancy(self,nrows=None):
        Nrows=self.shape[0]
        if not nrows: nrows=Nrows
        occ=numpy.zeros( (nrows,1) ) #use same shape as cm[0:nrows].sum(1)
        frame=0
        print 'frames processed...', ; sys.stdout.flush()
        for cm in self:
            if cm.data.any():
                occ+=cm[0:nrows].sum(1)
            frame+=1
            if not frame%1000: 
                print frame,'..', ; sys.stdout.flush()
        print '\n'
        return (occ/frame).reshape( (nrows,) ) #return as 1D array
    def ClusterRowsByKmeans(self,k,nrows=None,saturateOccupancy=False,jobName=None):
        """See Hoon10.pdf for cluster description
        For a given map NxP, add all the P column values to obtain a Nx1 vector.
        Each entry represents the number of contacts between a 'row' particle and all the 'column' particles.
        If we have M frames in the ContactMapList, we end up with a MxN matrix that we can cluster
        k: number of clusters/centroids
        nrows: cluster only the first nrows. If left to None, cluster all rows
        saturateOccupancy: if set to True, the Nx1 vectors will contain ony one or zeros
        jobName: save results of clusterin in a format suitable for Java TreeView"""
        from Pycluster import Record #http://pypi.python.org/pypi/Pycluster
        M = len(self)
        Nrows = self.shape[0]
        if not nrows: nrows = Nrows
        obs = numpy.zeros( (M,nrows) )
        iM = 0
        for cm in self:
            if cm.data.any(): 
                v = cm[0:nrows].sum(1)
                v = v.reshape( (nrows,) )
            else:
                v = numpy.zeros(nrows)
            if saturateOccupancy: v[v>0] = 1
            obs[iM] = v
            iM += 1
        record = Record()  #suitable to be viewed with Java Treeview
        record.data = obs
        record.uniqid = 'UNIQID'
        record.geneid = list(1+numpy.arange(M))
        record.expid = list(1+numpy.arange(nrows))
        clusterid,errors,nfound = record.kcluster(k,npass=25,method='a', dist='e')
        cdata, cmask = record.clustercentroids(clusterid=clusterid)
        if jobName: 
            record.save(jobName,clusterid,None)
        #find which frame is closest to each centroid
        nearestDistanceSquared = numpy.zeros(k)
        nearestMember = -numpy.ones(k,dtype='int')
        iM = 0
        for cm in self:
            if cm.data.any():
                v = cm[0:nrows].sum(1)
                v = v.reshape( (nrows,) )
            else:
                v = numpy.zeros(nrows)
            if saturateOccupancy: v[v>0] = 1
            icentroid = clusterid[iM]
            diff = v - cdata[ icentroid ]
            dd = numpy.square(diff).sum() #distance from v to the closest centroid
            if nearestMember[icentroid] < 0 or nearestDistanceSquared[icentroid] > dd:
                nearestMember[icentroid] = iM
                nearestDistanceSquared[icentroid] = dd
            iM += 1
        return clusterid,cdata,nearestMember

    def ResidenceTimes(self,nrows,FrameCutOff=0):
        """create a list of residence times
        FrameCutOff: append broken contacts only if bigger than FrameCutOff, or
            analogously, if time of broken contacts bigger then rtl.dt*FrameCutOff
        """
        rt=ResidenceTimesList(nrows)
        prev = self[0]       #first contact map
        R = copy(prev)             #store as initial residence times
        print 'frames processed...', ; sys.stdout.flush()
        iframe=1
        for curr in self[1:]:
            #positive elements of D indicate residence times of broken contacts
            D = R.multiply(prev-curr)
            R = R.multiply(curr) + curr    #update current residence times
            prev=curr                       
            rt.append(D,FrameCutOff) #append broken contacts only if bigger than FrameCutOff
            iframe += 1 
            if not iframe%1000: 
                print iframe,'..', ; sys.stdout.flush()
        print '\n'
        return rt
    def saveToHDF5(self,Hdf5Group):
        """save the list to a particular group in the HDF5 file"""
        idx=0
        for cm in self:
            subgroup = Hdf5Group.create_group('%05d'%idx)
            subgroup.attrs['index']=idx
            cm.saveToHDF5(subgroup)
            idx+=1
    def saveToFile(self,filename,mode='w',fmt='HDF5'):
        """save to file in the selected format"""
        if fmt=='HDF5':
            f = h5py.File(filename,mode)
            self.saveToHDF5(f['/'])
        else:
            sys.stderr.write('ERROR: format '+fmt+' not yet implemented\n')
    def loadFromFile(self,filename,fmt=None,offSet=0):
        """load list of contact maps from a file"""
        extDict={'hd5':'HDF5','txt':'DAT','dat':'DAT'}
        if not fmt:
            ext=filename.split('.')[-1]
            fmt=extDict[ext]
        if fmt=='HDF5':
            f = h5py.File(filename,'r')
            self.loadFromHDF5(f['/'],offSet=offSet)
            self.shape=self[0].shape
        else:
            print 'file format not understood'
    def loadFromHDF5(self,Hdf5Group,offSet=0):
        """load list of contact maps from a particular group in the HDF5 file"""
        groups=Hdf5Group.keys()
        groups.sort()
        for group in groups[offSet:]:
            self.append( load_csr_ContactMapFromHDF5(Hdf5Group[group]) )
        
            
    