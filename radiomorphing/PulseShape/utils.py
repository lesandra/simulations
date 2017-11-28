import numpy

def load_trace(directory, index, suffix=".trace"):
    """Load data from a trace file
    """
    path = "{:}/a{:}{:}".format(directory, index, suffix)
    with open(path, "r") as f:
        return numpy.array([map(float, line.split()) for line in f])

def getn(h):
    """Get the refractive index

       Reference:
        Zhaires (see email M. Tueros 25/11/2016)
    """
    # h in meters
    return 1. + 325E-06 * numpy.exp(-0.1218E-03 * h)

def getCerenkovAngle(h):
   """Get the Cerenkov angle
   """
   return numpy.arccos(1. / getn(h))
