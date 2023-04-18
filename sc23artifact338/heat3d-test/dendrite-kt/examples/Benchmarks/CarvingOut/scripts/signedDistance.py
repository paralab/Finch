import Multiprocessing
import multiprocessing
import trimesh
import sys
import numpy as np
import pickle
import warnings
from os import path
warnings.simplefilter("ignore", UserWarning)

def computeSignedDistance(mesh,trialNumber,numIter):
    fileName = "refineLevel"+str(trialNumber)+"/rank"+str(numIter)+".txt"
    if (not(path.exists(fileName))):
        return [-np.Inf,0,0]
    K = np.loadtxt(fileName)
    if((len(K.shape) == 1) and (K.shape[0] == 3)):
        K = K.reshape((1,3))

    LInfDistance = -np.Inf
    L2Distance = 0
    samplePointsPerIteration = 1000 # Compromise b/w memory and efficiency
    numIteration = int(np.ceil(K.shape[0]/samplePointsPerIteration))
    for i in range(0,numIteration):
        start:int = i*samplePointsPerIteration
        end: int = min((i+1)*samplePointsPerIteration,K.shape[0])
        distance = trimesh.proximity.signed_distance(mesh=mesh, points=K[start:end,:])
        LInfDistance = max(LInfDistance,np.max(distance))
        L2Distance = L2Distance + np.sum(distance**2)
    return [LInfDistance,L2Distance,K.shape[0]]



#
# # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    argc: int = len(sys.argv)

    if argc < 3:
        print('Usage: ', sys.argv[0],'stlFilename','RefinementIteration')
        exit(0)

    multiprocessing.set_start_method('spawn')
    stlFileName:str = sys.argv[1]
    refinementIteration = int(sys.argv[2])
    mesh = trimesh.load_mesh(stlFileName)

    totalNumProc: int = 192
    numberOfJobs:int = 16
    numIter:int = int(np.ceil(totalNumProc/numberOfJobs))
    # totalError = computeSignedDistance(mesh,0,80)
    # totalError = computeSignedDistance(mesh,0,18)
    # totalError = computeSignedDistance(mesh,0,19)
    LInfError = -np.Inf
    L2Error = 0
    numEntries = 0
    for iterID in range(0,numIter):
        mp = Multiprocessing.Multiprocessor()
        for proc in range(0, numberOfJobs):
            mp.run(computeSignedDistance,proc,mesh,refinementIteration,iterID*numberOfJobs+proc)
        error = mp.wait()
        for _err in error:
            if(LInfError < _err[0]):
                LInfError = _err[0]
            L2Error = L2Error + _err[1]
            numEntries = numEntries + _err[2]
        print('Numiter = ',iterID, 'of',numIter)

    print('RefinementIteration = ', refinementIteration, 'totalError',[LInfError,L2Error,numEntries])
    result = {'iter':refinementIteration,'LInfdistance':LInfError,'L2Distance':L2Error,'NumEntries':numEntries}
    dumpFileName = 'result'+str(refinementIteration)+'.dmp'
    pickle.dump(result,open(dumpFileName,'wb'))

