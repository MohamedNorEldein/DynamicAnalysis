import ctypes
import numpy as np
import matplotlib.pyplot as plt
from main import *



def calcResponceSpectrum(time, yg,Qu_COL, beta ,N):
    """
        beta is the ratio of the columns strength to the wall strength
    
    """

    sdof = SDOF(1, 0.0)   
    
    fwall = Force(1,2000)
    fcol = Force(1,Qu_COL)

    sdof.addForce(fwall)
    sdof.addForce(fcol)

    T0 = 0.01
    DT = 4.0

    sdof.set_Ground_Dis(  yg, time[1])
    sdof.set_Initial_Condition(0,0)
    T =[0] *N
    E =[0] *N
    R =[0] *N
    D =[0] *N

   
    for i in range(1,N):
        sdof.initBack()

        T[i]  = T0 +DT /N *(i)
        
        fwall.k =( (2*3.14/T[i])**2)*beta
        fcol.k = (1-beta)*( (2*3.14/T[i])**2)
        
        sdof.c = 0.0
        #print(f"dt = {time[1]} i = {i} c ={sdof.c} up = {f.up}")
        sdof.solve()
        R[i] = max(abs(sdof.get_a()[:-1]))
        D[i] = max(abs(sdof.get_U()[:-1]))
        
        E[i] = sdof.plasticEnergy()
        
        #draw(sdof.get_U(), sdof.get_Q(), "Q-u", "u","Q")
        #plt.show()
        
    return T,D, R, E


def calcEnergyeSpectrum_Bilinear(time, yg,ki, beta ,N):
    """
        beta is the ratio of the columns strength to the wall strength
    
    """

    sdof = SDOF(1, 0.0)   
    
    fwall = Force(ki*beta,2000)         # walls are responasable for final stiffness
    fcol = Force( ki*(1-beta),0)

    sdof.addForce(fwall)
    sdof.addForce(fcol)

    DQ = 30

    sdof.set_Ground_Dis(  yg, time[1])
    sdof.set_Initial_Condition(0,0)
    Q =[0] *N
    E =[0] *N
    R =[0] *N
    D =[0] *N

   
    for i in range(0,N):
        sdof.initBack()

        Q[i]  = DQ /N *(i)
        fcol.Qu = Q[i]
        
        sdof.solve()
        R[i] = max(abs(sdof.get_a()[:-1]))
        D[i] = max(abs(sdof.get_U()[:-1]))
        
        E[i] = sdof.plasticEnergy()
        
        '''draw(sdof.get_U(), sdof.get_Q(), "Q-u", "u","Q")
        plt.show()
        '''
    return Q,D, R, E


def maxARR(arr):
    imax=0
    for i in  range(0, len(arr)):
        if(arr[i]>arr[imax]):
            imax =i
    return imax


def minARR(arr):
    imin=0
    for i in  range(0, len(arr)):
        if(arr[i]<arr[imin]):
            imin =i
    return imin



if __name__ == "__main__":
    print("\n/*********************TEST***********************************/\n")
   
    t, sg = readData("SANTA MONICA.csv","t","a")
    tg, yg, vg,ag = change_time_step(t[1],sg,len(sg)*70,0,0)

    Rmin = []
    Qrmin = []
    Emax = []
    Qemax = []
    betas = []
    Kis = []
    slope = []

    for i in range(0, 1):
        for j in range(0,20):

            Ki = 1000
            beta = 0.05*j
            q1, D1,R1,E1 =  calcEnergyeSpectrum_Bilinear(tg,yg,Ki,beta,20)
            
            a = minARR(R1)
            Rmin.append( R1[a])
            Qrmin.append(q1[a])

            b = maxARR(R1)
            a = int((a+b)/2)
            s = (q1[a] - q1[a-1]) /(R1[a] - R1[a-1]) 
            slope.append(s)
            print(s)
            
            #draw(q1, R1,"Yermo","beta","s")
            #plt.show()
            a = maxARR(E1)
            Emax.append( E1[a])
            Qemax.append(q1[a])

            betas.append(beta)
            Kis.append(Ki)


            #print("-",end="", flush=True)
        print("\ni = ", i)


    draw(betas, slope,"Yermo","beta","s")
    plt.show()

    #data = np.array( [Rmin,Qrmin ,Emax ,Qemax,betas ,Kis ])
    #np.savetxt("csv_file.csv", data.transpose(), delimiter=",")

