# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 18:03:04 2016

@author: alib
"""

from earthenv import expatmos
from aerodynamicfunc import aerodynamicfunc
import math
from math import sqrt 

from numpy import*


class propulsion:
    def __init__(self):
        ''' Constructor for this class. '''
    
    def nominalPath(self,y,t,coast,vertAscent,pitchDuration,pitchRate,t_vertical,tVac,
                cl,cd,refA,massFlow, planet = 'earth'):
                
        density = []
        pressure = []
        time = []
        thrust = []
        lift = []
        AoA = []
        q = []    
        
        earthAtmosphere = expatmos()    # define environmental model 
        aero_forces = aerodynamicfunc() # create class to access aerodynamic methods
        if planet == 'earth':
            planetRadius = 6371000
            rotRate = 7.2722e-5
            gravity = 9.81
            slDensity = 1.2922  
            slPressure = 101325
            
        dydt = zeros_like(y) # create empty array dydt with equal length to u
        r = planetRadius + y[2]     # define orbital radius  
        
        # caluclate various aerodynamic values
        instantaneousPressure = earthAtmosphere.expdensity(y[2])       
        instantaneousDensity =  earthAtmosphere.expdensity(y[2])      
        instantaneousDynamicPressure = aero_forces.dynamicPressure(instantaneousDensity,y[3])
        instantaneousDrag = aero_forces.drag(instantaneousDynamicPressure,refA,cd)
        instantaneousLift = aero_forces.lift(instantaneousDynamicPressure,refA,cl)
        instantaneousMachNumber = aero_forces.machNumber(y[3],instantaneousPressure,instantaneousDensity) 
       
        # calculate wind velocity       
        w_e = 0
        w_n = 0
       
        # define rate of change of angle of attack during pitch over   
        if t > t_vertical and t <= t_vertical+pitchDuration:
            alpha = pi/2 - y[4] - pitchRate*(t - t_vertical)
        else:
            alpha = 0
            
        # set thrust = 0 if coast=True or set thrust to appropriate value given instantaneous pressure 
        if coast == False:
            T = self.thrustAtmosphere(tVac,instantaneousPressure,0.78)
        else:
            T = 0
            massFlow =  0
            
    
        
        # Equations of motion during vertical ascent - modified to avoid divide by
        # zero errors originiating with the original equations of motion    
        if vertAscent == True:
            
            alpha = 0   
            dydt[0] = 0              
            dydt[1] = 0
            dydt[2] = y[3]
            dydt[3] = 1/y[6] * (T-instantaneousDrag-y[6]*gravity) + r*rotRate**2
            dydt[4] = 0
            dydt[5] = 0
            dydt[6] = -massFlow
            
            
            # append instantaneous values to their respective list
            density.append(instantaneousDensity)
            pressure.append(instantaneousPressure)
            thrust.append(T)
            lift.append(instantaneousLift)
            AoA.append(alpha)
            time.append(t)
            q.append(instantaneousDynamicPressure)
            #print y[0],y[1],y[2],y[3],y[4],y[5],y[6]
        
            # print values of dydt and boolean conditions vertAscent and coast
         
         
        # assign appropriate thrust value given time of flight
        else:
            
        
            
            # powered ascent equations of motion without wind disturbance - note that wind will be added later
            dydt[0] = y[3]*cos(y[4])*cos(y[5])/(r*cos(y[1]))
            dydt[1] = y[3]*cos(y[4])*sin(y[5])/r
            dydt[2] = y[3]*sin(y[4])
            
            dydt[3] = 1/y[6] * (T*cos(alpha)-instantaneousDrag-y[6]*gravity*sin(y[4])+ \
                      w_e*cos(y[4])*cos(y[5])+w_n*cos(y[4])*sin(y[5]))+\
                      r*rotRate**2*cos(y[1])*(cos(y[1])*sin(y[4])-sin(y[1])*cos(y[4])*sin(y[5]))
                      
            dydt[4] = 1/(y[6]*y[3])*((T*sin(alpha)+instantaneousLift)- \
                      y[6]*gravity*cos(y[4])-w_n*sin(y[4])*sin(y[5])- \
                      w_e*sin(y[4])*cos(y[5]))+y[3]*cos(y[4])/r+2*rotRate*cos(y[1])*cos(y[5])+ \
                      r*rotRate**2/y[3]*cos(y[1])*(cos(y[1])*cos(y[4])+sin(y[1])*sin(y[4])*sin(y[5])) 
                      
            dydt[5] = -1/(y[6]*y[3]*cos(y[4]+.01))*((T*sin(alpha)+instantaneousLift)+ \
                      w_e*sin(y[5])-w_n*cos(y[5]))-y[3]/r*tan(y[1])*cos(y[4])*cos(y[5])+ \
                      2*rotRate*(cos(y[1])*tan(y[4])*sin(y[5])-sin(y[1]))- \
                      r*rotRate**2/(y[3]*cos(y[4]))*cos(y[1])*sin(y[1])*cos(y[5])   #0.0001 us there to avoid a divide by zero, prolly should go back and fix this at some point
            dydt[6] = -massFlow
            
           # print y[0],y[1],y[2],y[3],y[4],y[5],y[6]
            
        # append instantaneous values to their respective list
        density.append(instantaneousDensity)
        pressure.append(instantaneousPressure)
        thrust.append(T)
        lift.append(instantaneousLift)
        AoA.append(alpha)
        time.append(t)
        q.append(instantaneousDynamicPressure)       
     
        return dydt


    def thrustAtmosphere(self, vacuumThrust,atmosphericPressure,exitArea):
        return subtract(vacuumThrust,multiply(atmosphericPressure,exitArea))
        
    def burnTime(self, propellantMass, massFlow):
        
        if len(propellantMass) != len(massFlow):
            print ('Length of propellantMass and massFlow must be equal')
            quit   
            
        return (divide(propellantMass,massFlow))
            

        
        

        
    def massFlow(self, stageThrust, exVelocity, mixtureRatio):
        
        if len(stageThrust) != len(exVelocity) != len(mixtureRatio):
            print ('Length of stageThrust and exVelocity must be equal')
            quit
            
        
        propellantFlows = []
        fuelFlows = []
        oxidizerFlows = []
        
        for k in range(len(stageThrust)):
            
            massFlow = stageThrust[k]/exVelocity[k]
            oxidizerFlow = massFlow*mixtureRatio[k]/(mixtureRatio[k]+1)
            fuelFlow = massFlow/(mixtureRatio[k]+1)
           
            propellantFlows.append(massFlow)
            fuelFlows.append(fuelFlow)
            oxidizerFlows.append(oxidizerFlow)
            

        
        return propellantFlows,fuelFlows,oxidizerFlows
        

    def stageThrust(self, dryMass, maxG, planet = 'earth'):
        
        if len(dryMass) != len(maxG):
            quit
            print ('Length of dryMass and maxG must be equal')      
        
        if planet == 'earth':
            gravity = 9.8
        
        thrust = []
        for k in range(len(dryMass)):
           thrust.append(dryMass[k]*maxG[k]*gravity)
           
            
        return thrust
    
    def exRatio(self,gamma,pChamber,pExit):
        A = pow(2/(gamma+1),1/(gamma-1))
        B = pow(pChamber/pExit,1/gamma)
        C = (gamma+1)/(gamma-1)
        D = pow(pExit/pChamber,(gamma-1)/gamma)
        
        output = A*B/sqrt(C*(1-D))
        
        if output > 120: # prevent function from returning large expansion ratios
            return 120
        else:
            
            return output
        
    def thrustCoef(self,gamma,p1,p2,p3,expRatio,theoretical = False):
    
        A = 2*gamma**2/(gamma-1)
        B = pow(2/(gamma+1),(gamma+1)/(gamma-1))
        C = pow(p2/p1,(gamma-1)/gamma)
        D = (p2-p3)/p1*expRatio
        
        if theoretical == True :            
            cFactor = 1
        else:            
            cFactor = 0.97

        return  cFactor*math.sqrt(A*B*(1-C))+D
            
    def charVel(self,gamma,R,t1,theoretical = False):
        
        A = sqrt(gamma*R*t1)        
        B = math.sqrt(pow(2/(gamma+1),(gamma+1)/(gamma-1)))
        
        if theoretical == True :            
            cFactor = 1
        else:            
            cFactor = 0.95

        return A/(gamma*B)*cFactor
        
    def totalPressureRatio_inj2inlet(self,gamma,machInlet):
        A = 1+gamma*machInlet**2
        B = 1+(gamma-1)*machInlet**2*0.5
        C = gamma/(gamma-1)
        
        return A*pow(B,C)
        
    def staticPressureNozzleInlet(self,pInjector,gamma,machInlet):
        return pInjector/(1+gamma*machInlet**2)
        
    def staticPressureThroat(self, pChamper,gamma):
        return pChamper*pow(2/(gamma+1),(gamma/(gamma-1)))
        
    def areaRatio(self,gamma,pChamber,px):
        """ Function to determine the area ratio along a nozzle              
    
        Parameters
        ----------
        gamma : float
            ratio of specific heats [unitless]
        pChamber : float
            Stagnation pressure of combustion chamber [Pa]
        px : list, float
            Pressure along nozzle downwind from throat [Pa]        
            
        Returns
        -------
        y : {float}
            Area ratio [unitless] 
 
        See Also
        --------
        
        Notes
        -----
    
        
        Examples
        --------
        """  
        
        A = pow(2/(gamma+1),1/(gamma-1))
        B = pow(pChamber/px,1/gamma)
        C = (gamma+1)/(gamma-1)
        D = pow(px/pChamber,(gamma-1)/gamma)
        
        return A*B/sqrt(C*(1-D))
    
    def nozzleStagnationPressure(self,staticPressureNozzleInlet,machInlet,gamma):
        
        return pow(staticPressureNozzleInlet*(1+0.5*(gamma-1)*machInlet,gamma/gamma-1))
        
        
    def nozzlePressure(self,gamma,pChamber,areaRatio,tol = 1e-3):
        error = 100
        px = 10
        guess = 0
        tol = 1e-3
        while (error > tol):
            px += 1
            guess = self.areaRatio(gamma,pChamber,px)
            error = (guess-areaRatio)
        
        
        return px    
            
    
    def nozzlePressureProfile(self,nozzleExpansionRatio,gamma,pChamber):
        k = 1
        pressures = []
        expRatios = []
        while (k <= nozzleExpansionRatio):
            
            pressures.append(self.nozzlePressure(gamma,pChamber,k))
            expRatios.append(k)    
            k +=0.1
    
        return pressures,expRatios

    def inletTemperature(self,tChamber,gamma,machInlet):
        return tChamber/(1+0.5*(gamma-1)*machInlet**2)
        
    def nozzleTemperature(self,tChamber,px,pChamber,gamma):
        
        return multiply(tChamber,divide(px,pChamber)**((gamma-1)/gamma))
        
    def specificImpulse(self,cStar,cf):
        return cStar*cf/9.81
        
    def getExitRadius(self,throatRadius,expRatio):
     return math.sqrt(expRatio)*throatRadius
    
    def throatArea(self,thrust,cf,pChamber):
        ''' Description:  Calculates throat cross sectional area based on thrust, 
                      thrust coefficient and combustion chamber pressure         
        Units:
                thrust              [N]
                cf = thrust coefficient  [unitless]
                pChamber      [Pa]
        '''
        return thrust/(cf*pChamber)
    
    
    
    def area2Radius(self,a):
        return math.sqrt(a/math.pi)
    
    def frustrumConeVolume(self,h,R,r):
        return math.pi*h/3*(R**2+R*r+r**2)
        
    def chamberVol(self,r,R,l,fuel = 'RP1',ox = 'LOX'):
        
        if fuel == 'RP1' and ox == 'LOX':
            Lstar = 1.27
            
        throatArea = math.pi*R**2
        cor =  self.frustrumConeVolume(l,r,R)
        vol =  Lstar*throatArea
        cylVol = vol - cor
     
        return vol,cylVol
    
    def conicalLength(self,rt,ratio,r,a):
        ''' Description:  Determines the legnth ofa conical nozzle section
                  
        Units
        rt = throat radius              [m]
        r  = contour of a circular arc  [m]
        ratio = expansion or contraction ratio      [unitless]
        a  = half angle [radians]
        '''
        
        A = rt*(math.sqrt(ratio)-1)
        B = r*(1/math.cos(a)-1)
        return  (A+B)/math.tan(a)
        
    def chamberLength(self,volume,radius):
        return  volume/(math.pi*radius**2) 
    
    def exitVelocity(self,gamma,r,t1,p1,p2):
        A = divide(multiply(2,gamma),subtract(gamma,1))
        B = multiply(r,t1)
        C = subtract(1,pow(divide(p2,p1),divide(subtract(gamma,1),gamma)))
        
        return sqrt(multiply(multiply(A,B),C))
        
    def throatTemperature(self,chamberTemperature,throatPressure,chamberPressure,gamma):
        return pow(chamberTemperature*throatPressure/chamberPressure,(gamma-1)/gamma)
        
    def nozzleExitTemperature(self,chamberTemperature,exitPressure,chamberPressure,gamma):
        return pow(chamberTemperature*exitPressure/chamberPressure,(gamma-1)/gamma)