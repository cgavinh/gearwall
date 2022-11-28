import pickle

def pickle_save(dObj, sFilename):
    '''Given an object and a file name, write the object to the file using pickle.'''

    sFilename = sFilename + ".pickle"
    f = open(sFilename, "wb")
    p = pickle.Pickler(f)
    p.dump(dObj)
    f.close()

def pickle_load(sFilename):
    '''Given a file name, load and return the object stored in the file.'''

    sFilename = sFilename + ".pickle"
    f = open(sFilename, "rb")
    u = pickle.Unpickler(f)
    dObj = u.load()
    f.close()
    return dObj

def save(object, method, filename):
    if method == 'pickle':
        pickle_save(dObj=object, sFilename=filename)



def load(method, filename):
    if method == 'pickle':
        obj = pickle_load(sFilename=filename)
    if method == 'npz':
        obj = np.load(filename)
    return obj