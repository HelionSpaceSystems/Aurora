# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 18:03:04 2016

@author: alib
"""

from math import*
from numpy import*


class mass_estimates:
    def __init__(self):
        ''' Constructor for this class. '''
    
    
    def payloadFraction (self, mPayload, totalMass):
  
        return mPayload/totalMass
            
    def propellantMass(self, totalMass,emptyMass):
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
        
    def emptyMass(self,structuralRatio,totalMass):
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
            
    
    def optimumMassRatio(self, mPayload, exVelocity, lagrangeMultiplier, structuralRatio):        
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
        
    
    def burnoutVelocity(self, exVelocity, structuralRatio, payloadRatio):
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
            
            
            
            
    def circSpeed(self, orbitalAltitude, planet = 'earth'):
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
    
    def optimumStageMass(self, mPayload, orbitalAltitude, exVelocity, structuralRatio,
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
        vBurnout = self.burnoutVelocity(exVelocity, structuralRatio, payloadRatio)
    
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

    
    
        stageMass,massRatio = self.optimumMassRatio(mPayload,exVelocity,
                                               lagrangeMultiplier,structuralRatio)
        
        # check to find minimum stage mass is found
        for k in range(len(exVelocity)):
            check = lagrangeMultiplier*exVelocity[k]*(structuralRatio[k]*massRatio[k]-1)**2+\
                2*structuralRatio[k]*massRatio[k] - 1
            
            if check <= 0:
                print (' Unable to find minimum stage mass')
                quit            
        
        return stageMass,massRatio

    def growthFactor(self,mPaylado,stageMasses):
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