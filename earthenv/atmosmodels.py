# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 15:37:22 2015

@author: 
"""
from math import exp
    

class expatmos:
    def __init__(self):
        ''' Constructor for this class'''
        
    def expPressure(self,altitude):
        
        
        
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
            print "Altitude must be any real number grater than zero"
            
        
    
        
        
    def expdensity(self,altitude):
        
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
            print "Altitude must be any real number greater than zero"
        

       