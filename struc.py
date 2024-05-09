import numpy as np



class ShearBuilding:
    def __init__(self) -> None:
        self.m = []
        self.k = []
        self.M_mat =None
        self.K_mat = None
        self.eValues = None
        self.eVecs = None
        self.A = None
        self.B=None
        self.w = None
        self.T = None
    
    def add_floor(self, mass , stiffness):
        self.m.append(mass)
        self.k.append(stiffness)

    def build(self):
        n = len(self.m)
        self.K_mat = np.zeros((n,n))
        self.M_mat = np.zeros((n,n))

        for i in range(0, n-1):
            self.K_mat[i,i] = self.k[i]+self.k[i+1]
            self.K_mat[i,i+1] = -self.k[i+1]
            self.K_mat[i+1,i] = -self.k[i+1]
            self.M_mat[i,i] = self.m[i]

        self.K_mat[n-1,n-1] = self.k[n-1]
        self.M_mat[n-1,n-1] = self.m[n-1]
        
        mat = np.matmul( np.linalg.inv(self.M_mat), self.K_mat)
        self.eValues, self.eVecs = np.linalg.eig(mat)
        self.w = np.sqrt(self.eValues)
        self.T = 2 * np.pi / self.w


    def initialConditions(self, y0, u0):
        B = np.matmul( np.linalg.inv(self.eVecs),y0)
        A = np.matmul( np.linalg.inv(self.eVecs),u0/self.w)
        n = len(self.m)
        
        self.B = np.zeros((n,n))
        np.fill_diagonal(self.B,B)

        self.A = np.zeros((n,n))
        np.fill_diagonal(self.A,A)
        
    def calc_y(self,t):
        A = np.matmul(self.eVecs,self.A)
        B = np.matmul(self.eVecs,self.B)
        s = np.sin(self.w * t)
        c = np.cos(self.w *t)
        return  np.matmul(A,s) + np.matmul(B,c)

    def __str__(self):
        a =str( self.K_mat)
        a += "\n*********************************\n"
        a+= str(self.M_mat)
        a+= "\n Periodic times are:........ \n"
        a+= str(self.T)
        a+= "\neigen vectors are:........ \n"
        a+= str(self.eVecs)

        return a
        

if __name__ =='__main__':
                            
    building =  ShearBuilding()
     
    building.add_floor(3,2000)
    building.add_floor(2,1000)
    building.add_floor(4,3000)
    building.add_floor(1.5,1500)

    building.build()

    #print(building)
    building.initialConditions([0.1,0.15,0.2,0.3],[0,5,0,0])
    print(building.calc_y(5))