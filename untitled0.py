# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 09:36:49 2016

@author: alib
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 07:49:26 2016

@author: alib
"""
from prettytable import PrettyTable

from launch_vehicle import mass_estimates, propulsion
from earthenv import expatmos
from aerodynamicfunc import aerodynamicfunc
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from math import*
from numpy import*




stages = 3

earthAtmosphere = expatmos()    # define environmental model 
aero_forces = aerodynamicfunc() # create class to access aerodynamic methods
structure = mass_estimates()    # create class to access nass_estiame methods
propulsion = propulsion()

# payload input
mPayload = 151
payRatio = [0.2,0.2,0.28]


# propulsive incputs
mixtureRatio = [0, 2.3, 2.3]
exVelocity = [3.924,3.434,2.943]
gamma = [1.22,1.22,1.22] # lox/rp1 taken from liqued engine book page 85
R = [360, 360, 360] # lox/rp1 taken from liqued engine book page 85
tChamber = [3600,3600,3600] # lox/rp1 taken from liqued engine book page 85
tInj = tChamber # because we assume Minj = 0
pChamber = [7e6,7E6,7E6]
machInlet = [0.31,0.31,0.31] # page 6 liquid engine book typical value for gmama = 1.2
conHalfAngle = [40*math.pi/180,40*math.pi/180,40*math.pi/180] # typically 20 to 45 degrees
divHalfAngle = [15*math.pi/180,15*math.pi/180,15*math.pi/180] # typically 12 to 18 degrees
contractionRatio = [2.5,2.5,2.5]
nozzleDesignAltitude = [10E3,100E3,100E3] # stage 1, compromise high altitude and SL
                    

# orbital inputs
orbitalAltitude = 500

# structural inputs
structuralRatio = [0.2,0.2,0.2]
maxG = [1.0, 3, 6.] #g maximum acceptable acceleration per stage
refA = 5
engineCount = [1,1,9]


#aerodynamic inputs
cd = 0.2
cl = 0.2


# time events
t_vertical = 20
coastTimes = [10,10]
pitchDuration = 20
pitchRate = 0.001275
sep1Start2 = 5
sep2Start3 = 5
MECO2sep = 5
SECO2sep = 5


# launch site inputs
AvalonLong =    -53.10*pi/180
AvalonLat  = 46.64*pi/180

exVelocity = []
expRatio = []
charVel = []
totalPressureRatio_injector2Inlet = []
staticPressureNozzleIn = []
staticPressureThrt = []
throatArea = []
staticExitPressure = []
seaLevelThrustCoefficient = []
tInlet = []
idealThrustCoefficient = []
vacuumThrustCoefficient = []
isp = []
throatArea = []
throatRadius = []
throatTemp = []
exitRadius = []
staticPressureInjector = []
chamberRadius = []
conArc = []
divArc = []
conLength = []
divLength = []
chamberVolume = []
cylindricalSectionChamberVolume = []
cylChamberLength = []
totalChamberLength = []
exVelocity = []
thrusts = []
exitTemp = []


#for k in range(stages):
#    if k > 0 :
#        exVelocity.append(propulsion.exitVelocity(gamma[k],
#                                                  R[k],
#                                                  tChamber[k],
#                                                  pChamber[k],
#                                                  earthAtmosphere.expPressure(100)))
#                                                  
#    else:
#        exVelocity.append(propulsion.exitVelocity(gamma[k],
#                                                  R[k],
#                                                  tChamber[k],
#                                                  pChamber[k],
#                                                  earthAtmosphere.expPressure(0)))



exVelocity = [3.924,3.434,2.943]
stageMass,massRatio=structure.optimumStageMass(250, 600, exVelocity, 
                                               structuralRatio,payRatio)

emptyMasses = structure.emptyMass(structuralRatio,stageMass)   
propMass = structure.propellantMass(stageMass,emptyMasses)    


thrust = propulsion.stageThrust(emptyMasses,maxG)
thrustPerEngine = divide(thrust,engineCount)
propFlow,fuelFlow,oxidizerFlow = propulsion.massFlow(thrust,multiply(exVelocity,1000),mixtureRatio)

for k in range(stages):

    
    expRatio.append(propulsion.exRatio(gamma[k],pChamber[k],earthAtmosphere.expPressure(nozzleDesignAltitude[k])))
    charVel.append(propulsion.charVel(gamma[k],R[k],tChamber[k]))
    totalPressureRatio_injector2Inlet.append(propulsion.totalPressureRatio_inj2inlet(gamma[k],machInlet[k]))
    staticPressureInjector.append(multiply(pChamber[k],totalPressureRatio_injector2Inlet[k]))
    staticPressureNozzleIn.append(propulsion.staticPressureNozzleInlet(staticPressureInjector[k],gamma[k],machInlet[k]))
    staticExitPressure.append(propulsion.nozzlePressure(gamma[k],pChamber[k],expRatio[k]))    
    staticPressureThrt.append(propulsion.staticPressureThroat(pChamber[k],gamma[k]))
    vacuumThrustCoefficient.append(propulsion.thrustCoef(gamma[k],pChamber[k],staticExitPressure[k],earthAtmosphere.expPressure(100e3),expRatio[k]))
    
    tInlet.append(propulsion.inletTemperature(tChamber[k],gamma[k],machInlet[k]))
    

    if k > 0 :
        seaLevelThrustCoefficient.append('NA')
        isp.append(propulsion.specificImpulse(charVel[k],vacuumThrustCoefficient[k]))
        throatArea.append(propulsion.throatArea(thrustPerEngine[::-1][k],vacuumThrustCoefficient[k],pChamber[k]))


        
    else:
        seaLevelThrustCoefficient.append(propulsion.thrustCoef(gamma[k],pChamber[k],staticExitPressure[k],earthAtmosphere.expPressure(0),expRatio[k]))
        isp.append(propulsion.specificImpulse(charVel[k],seaLevelThrustCoefficient[k]))
        throatArea.append(propulsion.throatArea(thrustPerEngine[::-1][k],seaLevelThrustCoefficient[k],pChamber[k]))

        
    throatRadius.append(propulsion.area2Radius(throatArea[k]))
    exitRadius.append(propulsion.getExitRadius(throatRadius[k],15))
    chamberRadius.append(math.sqrt(contractionRatio[k])*throatRadius[k])
    conArc.append(throatRadius[k]*1.5)
    divArc.append(throatRadius[k]*0.328)
    conLength.append(propulsion.conicalLength(throatRadius[k],contractionRatio[k],conArc[k],conHalfAngle[k]))
    divLength.append(propulsion.conicalLength(throatRadius[k],contractionRatio[k],conArc[k],conHalfAngle[k]))
    Vc,cylVc = propulsion.chamberVol(chamberRadius[k],throatRadius[k],conLength[k])
    chamberVolume.append(Vc)
    cylindricalSectionChamberVolume.append(cylVc)
    cylChamberLength.append(propulsion.chamberLength(cylindricalSectionChamberVolume[k],chamberRadius[k]))
    totalChamberLength.append(conLength[k]+cylChamberLength[k])
    
    
    
    #throatTemp.append(propulsion.throatTemperature(tChamber[k],throatPressure[k],pChamber[k]))
    #nozzlePressure,expRatios = propulsion.nozzlePressureProfile(12,1.2,6e6)

    
    
    
payloadfrac = structure.payloadFraction(mPayload,(sum(stageMass)))
launchWetMass = sum(stageMass) + mPayload
burnTimes = propulsion.burnTime(propMass,propFlow)




structureTitles = ['Wet Mass [kg]',
                   'Dry Mass [kg]',
                   'Propellant Mass [kg]',
                   'Mass Ratio [kg]']


performanceTitles = ['Thrust [kN]',    
                     'Thrust/Engine [kN]',
                     'Sea Level Thrust Coefficient',
                     'Vacuum Thrust Coefficient',
                     'Charachteristic Velocity [m/s]',
                     'Exit Velocity [km/s]',
                     'Propellant Mass Flow [kg/s]',
                     'Oxidizer Mass Flow [kg/s]',
                     'Fuel Mass Flow [kg/s]',
                     'Burn Times [s]',
                     'Specific Impulse, Sea Level [s]']
                     
combustionChamberTitles = ['Chamber Radius [m]',
                           'Cylindrical Section Length [m]',
                           'Chamber Length [m]',
                           'Cylinder Section Volume [m^3]',
                           'Chamber Volume [m^3]',
                           'Throat Radius [m]',
                           'Throat Area [m^2]']
                           
nozzlePropertiesTitle = ['Expasion Ratio',
                         'Exit Radius [m]',
                         'Inlet Static Pressure [Pa]',
                         'Inlet Temperature [k]']     

injectorPropertiesTitle = ['Injector/Inlet Pressure Ratio',
                            'Injector Static Pressure [Pa]',
                            'Injector Temperature [K]']
                      
structureOutput = [stageMass[::-1],
                   emptyMasses[::-1],
                   propMass[::-1],
                   massRatio]
                   
combustionChamberOutput = [chamberRadius,
                           cylChamberLength,
                           totalChamberLength,
                           cylindricalSectionChamberVolume,
                           chamberVolume,
                           throatRadius,
                           throatArea]

performanceOutput = [divide(thrust[::-1],1000),
                    divide(thrustPerEngine[::-1],1000),
                    seaLevelThrustCoefficient,
                    vacuumThrustCoefficient,
                    charVel,
                    exVelocity,
                    propFlow[::-1],
                    oxidizerFlow[::-1],
                    fuelFlow[::-1],
                    burnTimes[::-1],
                    isp]

nozzlePropertiesOutput = [expRatio,
                          exitRadius,
                          staticPressureNozzleIn,
                          tInlet]

injectorPropertiesOutput = [totalPressureRatio_injector2Inlet,
                                staticPressureInjector,
                                tInj]
                            
performance = PrettyTable(["***Performance***",'Stage 1','Stage 2','Stage 3'])
performance.align[" "] = "l" # Left align city names

structure = PrettyTable(["***Structure***",'Stage 1','Stage 2','Stage 3'])
structure.align[" "] = "l" # Left align city names

combustionChamber = PrettyTable(['***Combustion Chamber***','Stage 1','Stage 2','Stage 3'])
combustionChamber.align[" "] = "l" # Left align city names

nozzleProperties = PrettyTable(['***Nozzle Properties***','Stage 1','Stage 2','Stage 3'])
nozzleProperties.align[" "] = "l" # Left align city names

injectorProperties = PrettyTable(['***Injector Properties***','Stage 1','Stage 2','Stage 3'])
injectorProperties.align[" "] = "l" # Left align city names


n = 0
for k in range(len(performanceTitles)):

    performance.add_row([performanceTitles[k],
                        performanceOutput[k][n],
                        performanceOutput[k][n+1],
                        performanceOutput[k][n+2]])


for k in range(len(structureTitles)):

    structure.add_row([structureTitles[k],
                      structureOutput[k][n],
                      structureOutput[k][n+1],
                      structureOutput[k][n+2]])

for k in range(len(combustionChamberTitles)):

    combustionChamber.add_row([combustionChamberTitles[k],
                      combustionChamberOutput[k][n],
                      combustionChamberOutput[k][n+1],
                      combustionChamberOutput[k][n+2]])

for k in range(len(nozzlePropertiesTitle)):

    nozzleProperties.add_row([nozzlePropertiesTitle[k],
                      nozzlePropertiesOutput[k][n],
                      nozzlePropertiesOutput[k][n+1],
                      nozzlePropertiesOutput[k][n+2]])


for k in range(len(injectorPropertiesTitle)):

    injectorProperties.add_row([injectorPropertiesTitle[k],
                      injectorPropertiesOutput[k][n],
                      injectorPropertiesOutput[k][n+1],
                      injectorPropertiesOutput[k][n+2]])         
                      
print 'Launch Wet Mass: ',launchWetMass,' kg'
print 'Payload Fraction: ',payloadfrac

print performance
print structure
print combustionChamber
print nozzleProperties
print injectorProperties





# Initialization of arrays to hold ouput from propulsion.nominalPath()

'''
density = []
pressure = []
time = []
lift = []
AoA = []
q = []
odeReturns = []

yLabel = ['Longitude [deg]','Lattitude [deg]','Altitude [km]',
                   'Velocity [km/s]','Flight Path Angle [deg]',
                   'Roll Angle [deg]','Mass [kg]','Density [kg/m^3]','','','','']

xLabel = ['Time','Time','Time','Time','Time','Time','Time','Time','Time','Time']






for k in range(len(stageMasses)*2+2):
    
    if k == 0:
        vertAscent = True
        coast = False
        yinit = [launchSiteCoords[0],launchSiteCoords[1],0,0,pi/2,0,launchWetMass]
        t0 = 0
        tf = t_vertical
        t = linspace(t0,tf,100)
        args = (coast,vertAscent,pitchDuration,pitchRate,t_vertical,thrust[-1],
                       cl,cd,refA,propFlow[-1])
        
        
        
    y = odeint(propulsion.nominalPath,yinit,t,args)




yLabel = ['Longitude [deg]','Lattitude [deg]','Altitude [km]',
                   'Velocity [km/s]','Flight Path Angle [deg]',
                   'Roll Angle [deg]','Mass [kg]','Density [kg/m^3]','','','','']
xLabel = ['Time','Time','Time','Time','Time','Time','Time','Time','Time','Time']    
fig, axs = plt.subplots(3,4, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = 0.5, wspace = 0.3)    
axs = axs.ravel()
for k in range(11):


    axs[k].set_ylabel(yLabel[k])
    axs[k].set_xlabel(xLabel[k])
    axs[k].plot(t,y[:,k])
'''

 

    


      