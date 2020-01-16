import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve as scipyFsolve
#get_ipython().run_line_magic('matplotlib', 'inline')
class field():
    def __init__(self, ncells, dx, nvars, bctype, additionalData = 1.4, bcs = np.array([0.0,0.0])):
        self.nvars = nvars
        self.nghost = 2 #Per side
        self.ncells = ncells
        self.dx = dx
        self.cc = np.zeros([ncells + 2*self.nghost,self.nvars])
        self.derivative = np.zeros([ncells + 2*self.nghost,self.nvars])
        self.cfVals_left = np.zeros([ncells+1,self.nvars])
        self.cfVals_right = np.zeros([ncells+1,self.nvars])
        self.flux_jphalf = np.zeros([ncells,self.nvars])
        self.flux_jmhalf = np.zeros([ncells,self.nvars])
        self.a = np.zeros(ncells+1)
        self.bc = bctype
        self.bcDirichlet = np.copy(bcs) #BCS must have the dimension 2*nvars
        self.gamma = additionalData
        
    def setInitCond(self, f, vals):
        for varidx in range(self.nvars):
            self.cc[self.nghost:-self.nghost,varidx] = f[varidx](vals)
            self.setGhostCellValues(varidx)
        self.computeCellFaceValues(self.dx)
        
    def computeDerivative(self, dx, theta = 1.0):
        self.derivative[1:-1] =  (1.0/dx)*vanLeer((self.cc[1:-1] - self.cc[0:-2]), (self.cc[2:] - self.cc[1:-1]))
    def setGhostCellValues(self, varidx):
        if (self.bc == 'periodic'):
            #print('Calling periodic BC')
            self.cc[:self.nghost,varidx] = self.cc[-2*self.nghost:-self.nghost,varidx]
            self.cc[-self.nghost:,varidx] = self.cc[self.nghost:2*self.nghost,varidx]
        elif (self.bc == 'dirichlet'):
            #print('Calling Dirichlet BC')
        #We will have problems if there are really strong gradients close to the boundary!
            self.cc[:self.nghost,varidx] = self.bcDirichlet[0,varidx]
            self.cc[-self.nghost:,varidx] = self.bcDirichlet[1,varidx]
        elif (self.bc == 'infiniteDomain'):
            #Has to be defined more mathematically
            raise ValueError('Under construction.....')
                
        else:
            raise ValueError('Unknown boundary condition')

    def computeCellFaceValues(self, dx):
        self.computeDerivative(dx)
        #print('rhs', self.cc[self.nghost-1:-self.nghost] + \
        #    0.5*dx*self.derivative[self.nghost-1:-self.nghost])
        self.cfVals_left[:] = self.cc[self.nghost-1:-self.nghost] + 0.5*dx*self.derivative[self.nghost-1:-self.nghost]
        self.cfVals_right[:] = self.cc[self.nghost:-self.nghost+1] - 0.5*dx*self.derivative[self.nghost:-self.nghost+1]
        
    def computeFluxes(self, dx, fluxFn):
        self.computeCellFaceValues(dx)
        
        self.a = self.computeLocalSpeed()
        
        for varidx in range(self.nvars):
            self.flux_jphalf[:,varidx] = 0.5*(fluxFn[varidx](self.cfVals_right[1:]) + fluxFn[varidx](self.cfVals_left[1:]) \
            					 - self.a[1:]*(self.cfVals_right[1:,varidx] -self.cfVals_left[1:,varidx]))
            
            self.flux_jmhalf[:,varidx] = 0.5*(fluxFn[varidx](self.cfVals_right[:-1]) + fluxFn[varidx](self.cfVals_left[:-1]) \
             					- self.a[:-1]*(self.cfVals_right[:-1,varidx] - self.cfVals_left[:-1,varidx]))
             					
             
        self.FinalFlux = -(self.flux_jphalf - self.flux_jmhalf)/dx
        
    def advance(self, dt):
        self.cc[self.nghost:-self.nghost] = self.cc[self.nghost:-self.nghost] + dt * self.FinalFlux
        for varidx in range(self.nvars):
            self.setGhostCellValues(varidx)
   
    def get_cc(self):
        return self.cc[self.nghost:-self.nghost]
    
    def printDeets(self):
        print(self.bc)
        
        
    def spectralRadius(self, cfVals):
        gasVel = cfVals[:,1]/cfVals[:,0]
        gasPressure = (self.gamma - 1.0)*(cfVals[:,2] - (cfVals[:,1]**2)/(2.0*cfVals[:,0]))
        soundVel = np.sqrt((self.gamma*gasPressure)/(cfVals[:,0]))
        return np.amax(np.abs([gasVel + soundVel, gasVel, gasVel - soundVel]),0)
        
    def computeLocalSpeed(self):
    #This is completely different from the previous case
        return np.amax(np.array([self.spectralRadius(self.cfVals_left), self.spectralRadius(self.cfVals_right)]),0)
    
    
    def returnPressure(self, uVals):
        return (self.gamma - 1.0)*(uVals[:,2] - (uVals[:,1]**2)/(2.0*uVals[:,0]))
    
    def returnVelocity(self, uVals):
        return uVals[:,1]/uVals[:,0]      

#------------------------------------------------------------------------------------------
class domain:
    def __init__(self, ncells, endpts):
        self.ncells = ncells
        self.xLeft = endpts[0]
        self.xRight = endpts[1]
        self.xMid = 0.5*(self.xLeft + self.xRight)
        self.dx = (self.xRight - self.xLeft)/self.ncells
        self.x_cc = np.linspace(self.xLeft+0.5*self.dx, self.xRight-0.5*self.dx, ncells)
#-------------------------------------------------------------------------------------------
                
#-------------------------------------------------------------------------------------------
def minmod(a, b):
    #This gives out a matrix of size (ncells*nvars)
    return 0.5 * (np.sign(a) + np.sign(b)) * np.amin([np.abs(a), np.abs(b)], axis=0)

def vanLeer(a,b):
    return (2.0*(np.maximum(0.0, a*b)))/(a + b + 1e-12)

def koren(a,b):
    aa = a*a
    ab = a*b
    phi2 = np.zeros(a.shape)
    flags = np.zeros(a.shape)
    flags[ab < 0] = 1.0
    flags[aa <= 0.25*ab] = 2.0
    flags[aa <= 2.5*ab] = 3.0
    phi2[flags == 1.0] = 0.0
    phi2[flags == 2.0] = a[flags == 2.0]
    phi2[flags == 3.0] = (1.0/6.0)*(b[flags == 3.0] + 2.0*a[flags == 3.0])
    phi2[flags == 0.0] = b[flags == 0.0]
    return phi2
#-------------------------------------------------------------------------------------------

def euler(Ncells, CFL, ics, tfinal, gamma):

    
    #Domain setup
    x = domain(Ncells, np.array([0.0,1.0]))
    gamma = 1.4
    
    
    def returnPressure(uVals):
        return (gamma - 1.0)*(uVals[:,2] - (uVals[:,1]**2)/(2.0*uVals[:,0]))
    
    def returnVelocity(uVals):
        return uVals[:,1]/uVals[:,0]

    def convertToConservatives(ics):
        #Input is [rhol,ul,pl],[rhor,ur,pr]
        c2l = ics[0,0]*ics[0,1];c2r = ics[1,0]*ics[1,1]
        c3l = ics[0,2]/(gamma-1.0) + 0.5*ics[0,0]*ics[0,1]**2
        c3r = ics[1,2]/(gamma-1.0) + 0.5*ics[1,0]*ics[1,1]**2
        return np.array([[ics[0,0],c2l,c3l],[ics[1,0],c2r,c3r]])


    def convertToPrimitives(ics):
        p2l = ics[0,1]/ics[0,0];p2r=ics[1,1]/ics[1,0]
        p3l = (gamma-1.0)*(ics[0,2] - 0.5*ics[0,1]**2/ics[0,0])
        p3r = (gamma-1.0)*(ics[1,2] - 0.5*ics[1,1]**2/ics[1,0])
        return np.array([[ics[0,0],p2l,p3l],[ics[1,0],p2r,p3r]])



    def exactSolution(ics, t):
        x0 = x.xMid
        rhoL = ics[0,0]; rhoR = ics[1,0]
        uL = ics[0,1]/rhoL; uR = ics[1,1]/rhoR
        pL = (gamma - 1.0)*(ics[0,2] - (ics[0,1]**2)/(2.0*ics[0,0]))
        pR = (gamma - 1.0)*(ics[1,2] - (ics[1,1]**2)/(2.0*ics[1,0]))
        aL = np.sqrt(gamma*pL/ics[0,0]); aR = np.sqrt(gamma*pR/ics[1,0])
        MachShock = lambda x: x - (1.0/x) - aL*((gamma+1)/(gamma-1))*(1.0 - ((pR/pL)*((2*gamma*x**2)/(gamma+1) - (gamma-1)/(gamma+1)))**((gamma-1)/(2*gamma)))

        Ms = scipyFsolve(MachShock, 6)
        #Computing the data in region 1
        p1 = pR*((2.0*gamma*Ms**2)/(gamma+1) - (gamma-1)/(gamma+1))
        u1 = (2/(gamma+1))*(Ms - 1.0/Ms)
        rho1 = rhoR/((2.0)/((gamma+1)*Ms**2) + (gamma-1)/(gamma+1))
        #Data in region 2
        p2 = p1; u2 = u1; rho2 = rhoL*(p2/pL)**(1.0/gamma); a2 = np.sqrt(gamma*(p2/rho2))
        #Data in region E
        x1 = x0 - aL*t; x2 = x0 + (u2 - a2)*t

        xELog = ((x1 < x.x_cc) & (x.x_cc < x2))
        uE = (2.0/(gamma+1))*(aL*xELog + (x.x_cc*(xELog) - x0*xELog)/t)
        aE = aL*xELog - (gamma-1)*uE*0.5
        pE = pL*(aE/aL)**((2.0*gamma)/(gamma-1))
        rhoE = (gamma*pE)/aE**2
        rhoE[np.isnan(rhoE)] = 0.0 #problems due to vectorization

        #Other positiions
        x3 = x0 + u2*t
        x4 = x0 + Ms*t #Position of shock(!)


        Pressure = (x.x_cc < x1)*pL + ((x1 < x.x_cc) & (x.x_cc < x2))*pE + ((x2 < x.x_cc) & (x.x_cc < x4))*p2 + (x.x_cc > x4)*pR 
        density = (x.x_cc < x1)*rhoL + ((x1 < x.x_cc) & (x.x_cc < x2))*rhoE + ((x2 < x.x_cc) & (x.x_cc < x3))*rho2 + ((x3 < x.x_cc) & (x.x_cc < x4))*rho1 + (x.x_cc > x4)*rhoR
        velocity = (x.x_cc < x1)*uL + ((x1 < x.x_cc) & (x.x_cc < x2))*uE + ((x2 < x.x_cc) & (x.x_cc < x4))*u2 + (x.x_cc > x4)*uR
        energy = (x.x_cc < x1)*uL + ((x1 < x.x_cc) & (x.x_cc < x2))*uE + ((x2 < x.x_cc) & (x.x_cc < x4))*u2 + (x.x_cc > x4)*uR
        xvals = np.copy(x.x_cc)
        #print('ShockLocation',x4)
        names = ['Density','Velocity','Pressure']
        #ax[0,0].plot(xvals, density,'r'), ax[0,0].set_xlabel('x'), ax[0,0].set_ylabel(names[0])
        #ax[0,1].plot(xvals, velocity,'b'), ax[0,1].set_xlabel('x'), ax[0,1].set_ylabel(names[1])
        #ax[1,0].plot(xvals, Pressure,'y'), ax[1,0].set_xlabel('x'), ax[1,0].set_ylabel(names[2])
    
    def makePlots(ufield, xvals, names, symbol):
        vals = ufield.get_cc() 
        ax[0,0].plot(xvals, vals[:,0],'r'+symbol), ax[0,0].set_xlabel('x'), ax[0,0].set_ylabel(names[0])
        ax[0,1].plot(xvals, returnVelocity(vals),'b'+symbol), ax[0,1].set_xlabel('x'), ax[0,1].set_ylabel(names[1])
        ax[1,0].plot(xvals, returnPressure(vals),'y'+symbol), ax[1,0].set_xlabel('x'), ax[1,0].set_ylabel(names[2])
    
    
    def f1(cfv):
        return cfv[:,1]
    
    def f2(cfv):
        return (cfv[:,1]**2/cfv[:,0]) + (gamma-1.0)*(cfv[:,2] - 0.5*cfv[:,1]**2/cfv[:,0])
    
    def f3(cfv):
        return (cfv[:,1]/cfv[:,0])*(cfv[:,2] + (gamma-1.0)*(cfv[:,2] - 0.5*cfv[:,1]**2/cfv[:,0]))
    
    fluxFn = [f1, f2, f3]
    nvars = 3
    #Initial conditions

    #ics = np.array([[1.0,0.0,2.5],[0.125,0.0,0.25]]) #Sod
    
    #ics = np.array([[0.445,0.311,8.928],[0.5,0.0,1.4275]]) #Lax
    
    #ics = convertToConservatives(np.array([[1.0,-2.0,0.8],[1.0,2.0,0.8]])) #123Problem - my code fails when -2.0
    
    #ics = convertToConservatives(np.array([[1.0,0.0,1000],[1.0,0.0,0.01]])) #Blast Wave- too slow/mesh not small
    
    #ics = convertToConservatives(np.array([[1.0,0.0,0.01],[1.0,0.0,100.0]])) #Woodward and Coallela
    
    ics = convertToConservatives(np.array([[5.9924,19.5975,460.894],[5.99242,-6.19633,46.095]])) #Two shocks colliding
    
    def ic1(c):
        return ics[0,0]*(c <= x.xMid) + (ics[1,0]*(c > x.xMid))
    
    def ic2(c):
        return ics[0,1]*(c <= x.xMid) + (ics[1,1]*(c > x.xMid))
    
    def ic3(c):
        return ics[0,2]*(c <= x.xMid) + (ics[1,2]*(c > x.xMid))
        
    initC = [ic1, ic2, ic3]
    
    u = field(x.ncells,x.dx, nvars, 'dirichlet', gamma, ics)
    u.setInitCond(initC, x.x_cc)
    fig, ax = plt.subplots(2,2, figsize=(12,12))
    makePlots(u, x.x_cc, ['Density','Velocity','Pressure'],'-')
    u.computeCellFaceValues(x.dx)
    t = 0.0
    tFinal = 0.035
    dfdu = u.computeLocalSpeed()
    dt = CFL * x.dx/np.amax(dfdu)
    n_steps = int(np.ceil(tFinal/dt))
    dt = tFinal / n_steps
    n_steps = int(np.ceil(tFinal/dt))
    for i in range(n_steps):
        old_cc = np.copy(u.cc)
        #Two step
        u.computeFluxes(x.dx, fluxFn)
        u.advance(dt)
        u.computeFluxes(x.dx, fluxFn)
        u.advance(dt)

        u.cc = 0.5*(old_cc + u.cc)
        t = t + dt
    exactSolution(ics, tFinal)
    makePlots(u, x.x_cc, ['Density','Velocity','Pressure'],'o')
    #exactSolution(ics, tFinal)
    plt.show()

