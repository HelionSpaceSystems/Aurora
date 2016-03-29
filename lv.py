# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 09:17:00 2016

@author: alib
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 18:03:04 2016

@author: alib
"""

from math import*
from numpy import*

import math
from math import sqrt 

from numpy import*

    
    
def payloadFraction (mPayload, totalMass):
  
    return mPayload/totalMass
        
def propellantMass(totalMass,emptyMass):
    """ Function to determine propellant mass of each stage in a launch vehicle
                
    Parameters
    ----------
    totalMass : array, float
        Array containing the total mass of each stage [kg]
        
    emptyMass : array, float
        Array containing the empty masses of each stage [kg]        
    Returns
    -------
    emptyMasses : {array, float}
        Empty mass of each stage [kg]

    
    See Also
    --------
    emptyMass: Function to determine empty mass of each stage in a launch vehicle        

    Notes
    -----
    Future planets to be added at a future date 
    
    Examples
    --------
    >>> propMass = propellantMass([475.8, 2217.1, 18753.5], [71.4, 443.4, 3750.7])
    >>> print (propMass)
    [71.4, 443.4, 3750.7] 
    """     
    
    propellantMass = []
    
    for k in range(len(totalMass)):
        propellantMass.append(totalMass[k] - emptyMass[k])
    
    return propellantMass
    
def emptyMass(structuralRatio,totalMass):
    """ Function to determine empty mass of each stage in a launch vehicle
                
    Parameters
    ----------
    structuralRatio : array, float
        Array containing structural ratio of each stage [unitless]
    totalMass : array, float
        Array containing the total mass [kg]
        
    Returns
    -------
    emptyMasses : {array, float}
        Empty mass of each stage [kg]

    
    See Also
    --------
    propellantMass: Function to determine propellant mass of each stage in a launch vehicle      

    Notes
    -----
    
    Examples
    --------
    >>> emptyMasses = emptyMass([0.15,0.2,0.2], 475.7, 2217.0, 18753.4)
    >>> print(emptyMasses)
    [71.36628000002011, 443.41822270378918, 3750.697085386098] 
    """     

    
    if len(structuralRatio) != len(totalMass):
        
        print
        ("The length of structuralRatio and totalMass must be equal")
        quit
    
    emptyMasses = []    
    for k in range(len(structuralRatio)):
        emptyMasses.append(structuralRatio[k]*totalMass[k])
    
    return emptyMasses
        

def optimumMassRatio(mPayload, exVelocity, lagrangeMultiplier, structuralRatio):        
    """ Returns optimum mass for each stage in a multi-stage launch vehicle              

    Parameters
    ----------
    mPayload : float
        Mass of payload [kg]
    exVelocity : array,float
        Array containing exhaust velocities of each stage, starting with stage one [km/s]
    lagrangeMultiplier : float
        Lagrange multiplier [unitless]
    structuralRatio : array, float
        Array containing structural ratio of each stage, starting with stage one [unitless]
    
    Returns
    -------
    stageMasses : {array, float}
        Array containing the mass of each stage, starting with the last stage [kg]
    stageRatios : {array. float}
        Array containing the stage ratios of each stage starting with the last stage [unitless]
    
    See Also
    --------
    optimumStageMass: Function to determine the optimum mass ratio of each stage in launch vehicle.
    
    Notes
    -----

    
    Examples
    --------
     
    """  
    
    
    # ensure exVelocity is the same length as structuralRatio
    if len(exVelocity) != len(structuralRatio):
        
        print ("The length of exVelocity and structuralRatio must be equal")
        quit

    # calculate mass ratio os each stage and append the result to stageRatios
    stageRatios = []
    stageMasses = []
    for x in range(len(exVelocity)):
        
        # divide equation into two intermediate variables, 'A' and 'B'
        A = exVelocity[x]*lagrangeMultiplier - 1
        B = exVelocity[x]*structuralRatio[x]*lagrangeMultiplier
        
        # append nth stage ratio to array stageRatios
        stageRatios.append(A/B)
    
   
    
    # calculate mass of nth stage and store in stageMasses
    stageRatios.reverse()   # reverse stageRatio
    structuralRatio.reverse()   # reverse structuralRatio
    stageMasses.append(mPayload)    #append mPayload to stageMasses 
    
    for x in range(len(exVelocity)):
        
        A = (stageRatios[x] - 1)*sum(stageMasses)
        B = 1 - stageRatios[x]*structuralRatio[x]
        stageMasses.append(A/B)

    stageMasses.pop(0)  #remove payload from returned result
    
    # re-arrange returns to begine with values for stage one

    return stageMasses,stageRatios
    

def burnoutVelocity(exVelocity, structuralRatio, payloadRatio):
    """     Function to calculate burnout velocity for multi stage launch vehicle              
    
    Parameters
    ----------
    exVelocity : array,float
        Array containing exhaust velocities of each stage, starting with stage one [km/s]
    structuralRatio : array, float
        Array containing structural ratio of each stage, starting with stage one [unitless]
    payloadRatio : array, float
        Array containing the payload ratio of each stage, starting with stage one [unitless]    
    Returns
    -------
    sum(burnoutVelocities) : {float}
        Sum of burnout velocitie of all stages [km/s]
    burnoutVelocities : {array. float}
        Array containing the burnout velocities of each stage starting with stage one [km/s]
    
    See Also
    --------
    circSpeed: Function to determine orbital speed of circular orbit

    
    Notes
    -----

    
    Examples
    --------
    >>> a = burnoutVelocity([3.924,3.434,2.943],[0.15,0.2,0.2],[0.15,0.2,0.2])
    >>> print(a)
    (11.545380378413567, [5.030310372902123, 3.5083503840248236, 3.0067196214866208])        
    """      
    
 
    
    if len(structuralRatio) != len(payloadRatio) != len(exVelocity):
        
        print ("The length of structuralRatio, exVelocity and payloadRatio must be equal")
        quit
    
    burnoutVelocities = []    
    for k in range(len(structuralRatio)):
        burnoutVelocities.append(-exVelocity[k]*math.log(structuralRatio[k]+\
        (1-structuralRatio[k])*payloadRatio[k],e))

    return sum(burnoutVelocities),burnoutVelocities
        
        
        
        
def circSpeed(orbitalAltitude, planet = 'earth'):
    """ Function to determine orbital speed of circular orbit
                
    Parameters
    ----------
    orbitalAltitude : float
        Altitude of desired orbit [km]
    planet : string
        Name of planet to be orbited     
    Returns
    -------
    y : {float}
        Orbital speed [km/s]

    
    See Also
    --------
    burnoutVelocity: Function to calculate burnout velocity for multi stage launch vehicle        

    Notes
    -----
    Future planets to be added at a future date 
    
    Examples
    --------
    >>> a = circSpeed(700, 'earth')
    >>> print(a)
    7.48331477355  
    """     

      
    earthRadius = 6371
    earhGravParameter = 398600
   
   
    if planet == 'earth':
        return sqrt(earhGravParameter/(earthRadius+orbitalAltitude))

def optimumStageMass(mPayload, orbitalAltitude, exVelocity, structuralRatio,
                     payloadRatio, tol = 10E-3):

    """ Function to determine the optimum mass ratio of each stage in launch vehicle.              

    Parameters
    ----------
    mPayload : float
        Mass of payload [kg]
    orbitalAltitude : float
        Altitude of desired orbit [km]        
    exVelocity : array,float
        Array containing exhaust velocities of each stage, starting with stage one [km/s]
    structuralRatio : array, float
        Array containing structural ratio of each stage, starting with stage one [unitless]
    payloadRatio : array, float
        Array containing the payload ratio of each stage, starting with stage one [unitless] 
        
    Returns
    -------
    stageMasses : {array, float}
        Array containing the mass of each stage, starting with the last stage [kg]
    stageRatios : {array. float}
        Array containing the stage ratios of each stage starting with stage one [unitless]
    
    See Also
    --------
    optimumMassRatio: Returns optimum mass for each stage in a multi-stage launch vehicle      
    
    Notes
    -----

    
    Examples
    --------
    >>> a = optimumStageMass(250, 1000, [3.924,3.434,2.943], [0.15,0.2,0.2], [0.15,0.2,0.38])
    >>> print(a)
    ([250, 475.77520000013413, 2217.0911135189458, 18753.485426930489], 
     [2.102751273746655, 2.5170055325091463, 3.7694179404133217])          
    """  
    
    if len(exVelocity) != len(structuralRatio) != len(payloadRatio):
        
        print ("The length of exVelocity and structuralRatio must be equal")
        quit
        
        
    # calculate burnout velocity for circular orbit
    vBurnout = burnoutVelocity(exVelocity, structuralRatio, payloadRatio)

    # iteratively solve for lagrange parameter
    lagrangeMultiplier = 0.7
    error = 1
    
    while(abs(error) >= tol):
        
        lagrangeMultiplier -= 0.0001
        
  
            
   
        A = sum([x*math.log(x*lagrangeMultiplier-1,e) for x in exVelocity])
        B = math.log(lagrangeMultiplier,e)*sum(exVelocity)      
        C = multiply(exVelocity,structuralRatio)
        D = [math.log(x,e) for x in C]
        E = sum(multiply(exVelocity,D))            
       
        guess = A - B - E   
        error = vBurnout[0] - guess



    stageMass,massRatio = optimumMassRatio(mPayload,exVelocity,
                                           lagrangeMultiplier,structuralRatio)
    
    # check to find minimum stage mass is found
    for k in range(len(exVelocity)):
        check = lagrangeMultiplier*exVelocity[k]*(structuralRatio[k]*massRatio[k]-1)**2+\
            2*structuralRatio[k]*massRatio[k] - 1
        
        if check <= 0:
            print (' Unable to find minimum stage mass')
            quit            
    
    return stageMass,massRatio

def growthFactor(mPaylado,stageMasses):
           ### INCOMPLETE ###
           
    """ Function to determine the optimum mass ratio of each stage in launch vehicle.              

    Parameters
    ----------
    mPayload : float
        Mass of payload [kg]
    launcWetMass : float
        Wet mass of launch vehicle [kg]
    stageMasses : list, float
        Mass of individual stages in the launch vehicle [km]        
        
    Returns
    -------
    y : {list, float}
        List containing the growth factor of each stage against the vehicle system weight 
 
    See Also
    --------
    
    Notes
    -----

    
    Examples
    --------
    >>> a = optimumStageMass(250, 1000, [3.924,3.434,2.943], [0.15,0.2,0.2], [0.15,0.2,0.38])
    >>> print(a)
    ([250, 475.77520000013413, 2217.0911135189458, 18753.485426930489], 
     [2.102751273746655, 2.5170055325091463, 3.7694179404133217])          
    """  
    
    y = [mPayload]+[stageWetMasses]
    launchWetMass = mPayload + sum(StageMasses)
    output = []
    
    for k in range(len(launchWetMass)):
        output.append(launchWetMass/(sum(stageMasses[k:]+mPayload)))
    return output


def nominalPath(y,t,coast,vertAscent,pitchDuration,pitchRate,t_vertical,tVac,
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


def thrustAtmosphere(vacuumThrust,atmosphericPressure,exitArea):
    return subtract(vacuumThrust,multiply(atmosphericPressure,exitArea))
    
def burnTime(propellantMass, massFlow):
    
    if len(propellantMass) != len(massFlow):
        print ('Length of propellantMass and massFlow must be equal')
        quit   
        
    return (divide(propellantMass,massFlow))
        

    
    

    
def massFlow(stageThrust, exVelocity, mixtureRatio):
    
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
    

def stageThrust(dryMass, maxG, planet = 'earth'):
    
    if len(dryMass) != len(maxG):
        quit
        print ('Length of dryMass and maxG must be equal')      
    
    if planet == 'earth':
        gravity = 9.8
    
    thrust = []
    for k in range(len(dryMass)):
       thrust.append(dryMass[k]*maxG[k]*gravity)
       
        
    return thrust

def exRatio(gamma,pChamber,pExit):
    A = pow(2/(gamma+1),1/(gamma-1))
    B = pow(pChamber/pExit,1/gamma)
    C = (gamma+1)/(gamma-1)
    D = pow(pExit/pChamber,(gamma-1)/gamma)
    
    output = A*B/sqrt(C*(1-D))
    
    if output > 120: # prevent function from returning large expansion ratios
        return 120
    else:
        
        return output
    
def thrustCoef(gamma,p1,p2,p3,expRatio,theoretical = False):

    A = 2*gamma**2/(gamma-1)
    B = pow(2/(gamma+1),(gamma+1)/(gamma-1))
    C = pow(p2/p1,(gamma-1)/gamma)
    D = (p2-p3)/p1*expRatio
    
    if theoretical == True :            
        cFactor = 1
    else:            
        cFactor = 0.97

    return  cFactor*math.sqrt(A*B*(1-C))+D
        
def charVel(gamma,R,t1,theoretical = False):
    
    A = sqrt(gamma*R*t1)        
    B = math.sqrt(pow(2/(gamma+1),(gamma+1)/(gamma-1)))
    
    if theoretical == True :            
        cFactor = 1
    else:            
        cFactor = 0.95

    return A/(gamma*B)*cFactor
    
def totalPressureRatio_inj2inlet(gamma,machInlet):
    A = 1+gamma*machInlet**2
    B = 1+(gamma-1)*machInlet**2*0.5
    C = gamma/(gamma-1)
    
    return A*pow(B,C)
    
def staticPressureNozzleInlet(pInjector,gamma,machInlet):
    return pInjector/(1+gamma*machInlet**2)
    
def staticPressureThroat( pChamper,gamma):
    return pChamper*pow(2/(gamma+1),(gamma/(gamma-1)))
    
def areaRatio(gamma,pChamber,px):
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

def nozzleStagnationPressure(staticPressureNozzleInlet,machInlet,gamma):
    
    return pow(staticPressureNozzleInlet*(1+0.5*(gamma-1)*machInlet,gamma/gamma-1))
    
    
def nozzlePressure(gamma,pChamber,ratio,tol = 1e-3):
    error = 100
    px = 10
    guess = 0
    tol = 1e-3
    while (error > tol):
        px += 1
        guess = areaRatio(gamma,pChamber,px)
        error = (guess-ratio)
    
    
    return px    
        

def nozzlePressureProfile(nozzleExpansionRatio,gamma,pChamber):
    k = 1
    pressures = []
    expRatios = []
    while (k <= nozzleExpansionRatio):
        
        pressures.append(nozzlePressure(gamma,pChamber,k))
        expRatios.append(k)    
        k +=0.1

    return pressures,expRatios

def inletTemperature(tChamber,gamma,machInlet):
    return tChamber/(1+0.5*(gamma-1)*machInlet**2)
    
def nozzleTemperature(tChamber,px,pChamber,gamma):
    
    return multiply(tChamber,divide(px,pChamber)**((gamma-1)/gamma))
    
def specificImpulse(cStar,cf):
    return cStar*cf/9.81
    
def getExitRadius(throatRadius,expRatio):
 return math.sqrt(expRatio)*throatRadius

def throatArea(thrust,cf,pChamber):
    ''' Description:  Calculates throat cross sectional area based on thrust, 
                  thrust coefficient and combustion chamber pressure         
    Units:
            thrust              [N]
            cf = thrust coefficient  [unitless]
            pChamber      [Pa]
    '''
    return thrust/(cf*pChamber)



def area2Radius(a):
    return math.sqrt(a/math.pi)

def frustrumConeVolume(h,R,r):
    return math.pi*h/3*(R**2+R*r+r**2)
    
def chamberVol(r,R,l,fuel = 'RP1',ox = 'LOX'):
    
    if fuel == 'RP1' and ox == 'LOX':
        Lstar = 1.27
        
    throatArea = math.pi*R**2
    cor =  frustrumConeVolume(l,r,R)
    vol =  Lstar*throatArea
    cylVol = vol - cor
 
    return vol,cylVol

def conicalLength(rt,ratio,r,a):
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
    
def chamberLength(volume,radius):
    return  volume/(math.pi*radius**2) 

def exitVelocity(gamma,r,t1,p1,p2):
    A = divide(multiply(2,gamma),subtract(gamma,1))
    B = multiply(r,t1)
    C = subtract(1,pow(divide(p2,p1),divide(subtract(gamma,1),gamma)))
    
    return sqrt(multiply(multiply(A,B),C))
    
def throatTemperature(chamberTemperature,throatPressure,chamberPressure,gamma):
    return pow(chamberTemperature*throatPressure/chamberPressure,(gamma-1)/gamma)
    
def nozzleExitTemperature(chamberTemperature,exitPressure,chamberPressure,gamma):
    return pow(chamberTemperature*exitPressure/chamberPressure,(gamma-1)/gamma)
    

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 18:06:36 2015

@author: alib
"""



def machNumber(v,p,rho):
    
    try:
        v = float(v)
        p = float(p)
        rho = float(rho)
        
        v >= 0
        p >= 0
        rho >= 0
        
        return v/sqrt(1.4*p/rho)
        
    except ValueError:
        print ("All inputs must be real numbers greater than zero")
        

def dynamicPressure(rho,v):
    
    try:
        rho = float(rho)
        v = float(v)

        
        rho >= 0
        v >= 0

        
        return 0.5*rho*v**2
        
    except ValueError:
        print ("All unputs must be real numbers greater than zero")
        

def drag(q,s,cd):
    
    inputs = [q,s,cd]
    
    
    for x in inputs:
        try:
            x = float(x)
            x >= 0
        except ValueError:
            print ("All inputs must be real numbers greater than or equal to zero")
            
    return q*s*cd

def lift(q,s,cl):
    
    inputs = [q,s,cl]
    
    
    for x in inputs:
        try:
            x = float(x)
            x >= 0
        except ValueError:
            print ("All inputs must be real numbers greater than or equal to zero")
            
    return q*s*cl




        
            
def expPressure(altitude):   
    # define constants
    seaLevelPressure = 101326   # sea level pressure [Pa]
    scaleHeight = 8240          # scale height [m]
    
    # validate user input
    try:
        
        altitude = float(altitude) # reject any strings
        altitude >= 0              # reject nonsensical altitudes
        
        #calculate pressure given altitude
        pressure = seaLevelPressure*exp(-altitude/scaleHeight)
        return pressure
        
    except ValueError:
        print ("Altitude must be any real number grater than zero")
        
      
    
def expdensity(altitude):
    
    # define constants    
    seaLevelDensity = 1.23  # sea level density kg/m^3
    scaleHeight = 8240      # scale height [m]
    
    # validate user inpit 
    try:
        
        altitude = float(altitude)
        altitude >=0 
        
        # calculate density given altitude
        density = seaLevelDensity*exp(-altitude/scaleHeight)
        return density 

        
    except ValueError:
        print ("Altitude must be any real number greater than zero")
    



def coe2sv(coe):
    """ Function to determine the area ratio along a nozzle              

    Parameters
    ----------
    mu : gravitational parameter (km^3/s^2)
    coe : list of orbital elements [h,w,RA,incl,w,TA]
        h: angular momentum (km^2/s)
        e: eccentricity
        RA: right ascension of the ascending node (rad)
        incl: inclination (rad)
        w: argument of perigee (rad)
        TA: true anomaly (rad)
    R3_w: rotation matrix abour the z axis through the angle w
    R1_i: rotationa matrix about the x axis through the angle i
    R3_W: rotation matrixabout the z axis through the angle RA
    Q_pX: matric of the transformation from perifocal to geocentric equatorial frame
    rp: position vector in the perifocal frame (km)
    vp: velocity vector in the perifocal frame (km/s)

        
    Returns
    -------
    sv : {list} [r,v]
        r: position vector in the geocentric equatorial frame (km)
        v - velocity vector in the geocentric equatorial frame (km/s)
 
 
    
    See Also
    --------
    
    Notes
    -----

    
    Examples
    --------
    """  
    
    R3_w = np.matrix('cos(RA) sin(RA) 0;-sin(RA) cos(RA) 0;0 0 1')
    print (R3_w)
    