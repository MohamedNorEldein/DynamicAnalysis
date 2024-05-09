import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import ctypes
import numpy as np
import scipy

def draw(x_values,y_values ,label1 = '', x = '',  y= ''):
        # Plot the data
    plt.plot(x_values, y_values, label=label1)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.legend()  # Display legend if multiple plots are present
    plt.grid(True)

    return


def readData(txt, a, b):
        # Lists to store x and y values
    x_values = []
    y_values = []

    # Read data from CSV file
    print("reading CSV file ")
    with open(txt, 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        print(csv_reader.fieldnames)
        for row in csv_reader:
            x_values.append(float(row[a]))
            y_values.append(float(row[b]))
    
    print("finished reading CSV file ")

    return x_values,y_values


def generate(count ,Amplitude, w,dt):
    yg = np.zeros(count, dtype=np.float64)
    vg = np.zeros(count, dtype=np.float64)
    ag = np.zeros(count, dtype=np.float64)

    t = np.zeros(count, dtype=np.float64)
    for i in range(0,count):
        t[i] = dt*i
        yg[i] = Amplitude * np.sin(w * i *dt)/w/w
        vg[i] = Amplitude * np.cos(w * i *dt)/w
        ag[i] = -Amplitude * np.sin(w * i *dt)


    return t,yg,vg,ag


#--------------------------------------
# Load the DLL
dll_path = ".\\build\\Release\\MyProject.dll"  # Replace with the actual path to your DLL
lib = ctypes.CDLL(dll_path)


class Force(ctypes.Structure):
    _fields_ = [("k", ctypes.c_double),
                ("Qu", ctypes.c_double),
                ("Ep", ctypes.c_double),
                ("up", ctypes.c_double)
                ]

    def __init__(self, k, Qu, Ep=0.0, up =0.0):
        self.k = k
        self.Qu = Qu
        self.Ep = Ep
        self.up = up

    def __call__(self, y) :
        return lib.forceFunc(self,ctypes.c_double(y))

   

class SDOF(ctypes.Structure):
    _fields_ = [
        ("yg", ctypes.POINTER(ctypes.c_double)),
        ("a", ctypes.POINTER(ctypes.c_double)),
        ("v", ctypes.POINTER(ctypes.c_double)),
        ("y", ctypes.POINTER(ctypes.c_double)),
        ("dt", ctypes.c_double),
        ("C", ctypes.c_double),
        ("m", ctypes.c_double),
        ("Q", ctypes.POINTER(ctypes.c_double)),
        ("u", ctypes.POINTER(ctypes.c_double)),
        ("forces", ctypes.POINTER(Force) * 10),
        ("fCount", ctypes.c_size_t),
        ("count", ctypes.c_size_t)
       
    ]

    def __init__(self, mass, damping):
        super().__init__()
        self.ag = None
        self.a = None
        self.v = None
        self.y = None
        
        self.Q = None
        self.u = None

        self.dt = 0
        self.Ep = 0
        self.c = damping
        self.m = mass
        self.fCount = 0
        self.count = 0

    def addForce(self ,force:Force):
        lib.addForce(ctypes.byref(self), ctypes.byref(force))

    def calcForce(self,  dx):
        return lib.calcForceSDOF(ctypes.byref(self), dx)

    def set_Ground_Dis(self, ag, dt):
        ag_ptr = np.ctypeslib.as_ctypes(ag)
        lib.setGroundDis(ctypes.byref(self), ag_ptr, len(ag), dt)

    def set_Initial_Condition(self, y0, v0):
        lib.setInitialCondition(ctypes.byref(self),y0,v0)

    def solve(self):
        lib.solve(ctypes.byref(self))

    def plasticEnergy(self):
        Ep =0
        for i in range(0,self.fCount):
            Ep += self.forces[i].contents.Ep
            
        return Ep

    def initBack(self):
        for i in range(0,self.fCount):
            self.forces[i].contents.up = 0
            self.forces[i].contents.Ep = 0
            

    def get_a(self):
        return  np.ctypeslib.as_array( self.a,(self.count,))

    def get_v(self):
        return  np.ctypeslib.as_array( self.v,(self.count,))
        
    def get_y(self):
       return  np.ctypeslib.as_array( self.y,(self.count,))

    def get_Q(self):
        return  np.ctypeslib.as_array( self.Q,(self.count,))

    def get_U(self):
        return  np.ctypeslib.as_array( self.u,(self.count,))
        

# forceFunc
lib.forceFunc.restype = ctypes.c_double
lib.forceFunc.argtypes = [ctypes.POINTER(Force), ctypes.c_double]

# createSDOF
lib.createSDOF.restype = ctypes.POINTER(SDOF)
lib.createSDOF.argtypes = [ctypes.c_double, ctypes.c_double]

# addForce
lib.addForce.restype = ctypes.c_int
lib.addForce.argtypes = [ctypes.POINTER(SDOF), ctypes.POINTER(Force)]

# calcForceSDOF
lib.calcForceSDOF.restype = ctypes.c_double
lib.calcForceSDOF.argtypes = [ctypes.POINTER(SDOF), ctypes.c_double, ctypes.c_double]

# setGroundDis
lib.setGroundDis.restype = None
lib.setGroundDis.argtypes = [ctypes.POINTER(SDOF), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
                                      ctypes.c_double]

#CAPI void setInitialCondition(SDOF *sys,double y0,double v0)
lib.setInitialCondition.restype = None
lib.setInitialCondition.argtypes =[ctypes.POINTER(SDOF), ctypes.c_double,ctypes.c_double]


# solve
lib.solve.restype = None
lib.solve.argtypes = [ctypes.POINTER(SDOF)]

# Define utility functions
def createSDOF(m, c):
    sdof_ptr = lib.createSDOF(m, c)
    return sdof_ptr.contents

def get_data(s, count):
    if s is not None:
        return [s[i] for i in range(count)]
    else:
        return []


# Define the function signature for the C function
func_changeTimeStep = lib.changeTimeStep
func_changeTimeStep.argtypes = [ctypes.c_double, np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_size_t,
                                np.ctypeslib.ndpointer(dtype=np.float64),np.ctypeslib.ndpointer(dtype=np.float64),
                                np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_size_t]
func_changeTimeStep.restype = ctypes.c_double

# Define a Python wrapper for the C function
def change_time_step(dt_old, s,new_count, y0, v0):
    count = len(s)
    print(new_count)
    y = np.zeros(new_count, dtype=np.float64)
    v = np.zeros(new_count, dtype=np.float64)
    a = np.zeros(new_count, dtype=np.float64)

    y[0] = y0
    v[0] = v0

    # Call the C function
    result = func_changeTimeStep(ctypes.c_double(dt_old), np.array(s,dtype=np.float64), count, a, v, y, ctypes.c_size_t(new_count))
    
    t = np.zeros(new_count)

    for i in range(0, new_count):
        t[i] = i*result

    print(count , new_count)
    return t, y, v, a
