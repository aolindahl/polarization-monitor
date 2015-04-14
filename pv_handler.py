import pyca
from Pv import Pv
from numpy import nan as npNan

pvDefinitions = [['intensities', 'AMO:CB:DET', 2040],
                ['spectrum', 'AMO:CB:ES', 1026],
                ['energy', 'AMO:CB:PE', 360],
                ['polarization', 'AMO:CB:POL', 600],
                ['ebeam', 'AMO:CB:EBEAM', 240]]


class PvHandler:
    def __init__(self, timeout=1.0):
        # The definition of the pv names
        self._pvNames = {}
        self._pvSize = {}
        for pvDef in pvDefinitions:
            self._pvNames[pvDef[0]] = pvDef[1]
            self._pvSize[pvDef[0]] = pvDef[2]



        # Make the PV objects
        self._pv = {}
        for k, v in self._pvNames.iteritems():
            self._pv[k] = Pv(v)


        # Open connection to the PVs 
        try:
            for pv in self._pv.itervalues():
                pv.connect(timeout=timeout)
        except pyca.pyexc, e:
            print 'ERROR: Failed to connect:', e
            raise
        except pyca.caexc, e:
            print 'ERROR: Channel access error:', e
            raise

        self._timeout = timeout


    def __del__(self):
        for pv in self._pv.itervalues():
            pv.disconnect()


    def assignData(self, intensities=None, spectrum=None, energy=None,
            polarization=None, ebeam=None, verbose=False):
        input = {'intensities': intensities,
                'spectrum' : spectrum,
                'energy' : energy,
                'polarization': polarization,
                'ebeam': ebeam}
        
        for k, data in input.iteritems():
            if data is not None:
                if len(data) < self._pvSize[k]:
                    tupleData = tuple(data) + (npNan,)*(self._pvSize[k] - len(data))
                elif len(data) > self._pvSize[k]:
                    tupleData = tuple( data[-self._pvSize[k]:] )
                else:
                    tupleData = tuple( data )

                try:
                    if verbose:
                        print 'tupleData have type {} and size {}.'.format(
                                type(tupleData), len(tupleData))
                        print 'First objects with types:',
                        for i in range(5):
                            print '\t{}:{}'.format(tupleData[i],
                                    type(tupleData[i]))
                        print 'tupleData: {}'.format(tupleData)
                    self._pv[k].put( value=tupleData, timeout=self._timeout )
                except:
                    print 'ERROR: Failed to put data into pv {}.'.format(
                            self._pvNames[k])
                    raise

    def flushData(self):
        pyca.flush_io()



if __name__ == '__main__':
    import numpy as np
    import time

    pvHandler = PvHandler(timeout=1.0)
    
    
    t = time.time()
    for i in range(10) :
        print i
    
        data = {}
        for pvDef in pvDefinitions:
           data[pvDef[0]] = np.random.rand(pvDef[2])

        pvHandler.assignData(**data)
        pvHandler.flushData()
        print 'Processing time is {} s.'.format(time.time()-t)
        time.sleep(1.0)
        t = time.time() 



    print "Done!"
