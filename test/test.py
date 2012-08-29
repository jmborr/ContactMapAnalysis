import argparse
import sys,os
import numpy
from pdb import set_trace as trace

import ContactMapAnalysis as CMA
import ContactMapAnalysisAPI as CMAPI

def ContactMap(*kargs,**kwargs):
    
    jobList = ['generate residue ContactMaps',
               'number of atomic contacts',
               'occupancy',
               ]
    JOB = kwargs['job']  #what job to do
    if JOB == 'generate residue ContactMaps':
        """calculate the protein-protein residue-based contact maps"""
        #variables defining the ContactMapProtocol
        rowSelection = kwargs['rowSelection'] #protein residues
        colSelection = rowSelection           #self-contact map
        cutOff=float( kwargs['cutOff'] )      #atomic cutOff
        pargs = {'rowSelection':rowSelection,'colSelection':colSelection,'cutOff':cutOff,'byres':True}
        protocol = CMAPI.ContactMapProtocol(**pargs)
        #generate the ContactMapList with the trajectory and protocol
        PSFile = kwargs['PSFile']      #topology
        trajfile = kwargs['trajfile']  #trajectory
        cml = CMAPI.GenerateContactMapList(PSFile,trajfile,protocol) #contactMapList
        cml.saveToFile(kwargs['outf'],fmt='HDF5') #save to file in HDF5 format 
    elif JOB == 'number of atomic contacts':
        """calculate number of protein-protein contacts for each frame"""
        cml = CMA.ContactMapList()
        cml.loadFromFile(kwargs['inFile'])
        outStr = "#frame #contacts\n"
        frame = 1
        for cm in cml:
            outStr += '%5d %3d\n'%(frame,cm.nnz)
            frame += 1
        open(kwargs['outFile'],'w').write(outStr)
    else:
        print 'Job not recognized. Nothing to do'



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='service provider for dhfr_solv project')
    parser.add_argument('service',help='requested service, the name of a function defined in this module')
    parser.add_argument('--kargs',help='required arguments of service. Ex: "arg1,arg2,arg3"')
    parser.add_argument('--kwargs',help='optional arguments of service. Ex: "arg1=val1,arg2=val2"')
    args = parser.parse_args()
    service=args.service
    reqargs=[]
    if args.kargs: reqargs=args.kargs.split(',')
    optargs={}
    trace()
    if args.kwargs: optargs=dict( map( lambda x: x.split('='), args.kwargs.split(',') ) )
    exitCode=0
    try:
        locals()[service](*reqargs,**optargs)
    except:
        sys.stderr.write('Error: service '+service+' not found or error in service\n')
        exitCode=1
    sys.exit(exitCode)