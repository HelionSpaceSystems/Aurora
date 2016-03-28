# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 18:06:36 2015

@author: alib
"""

from math import sqrt

class aerodynamicfunc:
    def __init__(self):
        ''' Constructor for this class'''
    
    def machNumber(self,v,p,rho):
        
        try:
            v = float(v)
            p = float(p)
            rho = float(rho)
            
            v >= 0
            p >= 0
            rho >= 0
            
            return v/sqrt(1.4*p/rho)
            
        except ValueError:
            print "All inputs must be real numbers greater than zero"
            

    def dynamicPressure(self,rho,v):
        
        try:
            rho = float(rho)
            v = float(v)

            
            rho >= 0
            v >= 0

            
            return 0.5*rho*v**2
            
        except ValueError:
            print "All unputs must be real numbers greater than zero"
            
    
    def drag(self,q,s,cd):
        
        inputs = [q,s,cd]
        
        
        for x in inputs:
            try:
                x = float(x)
                x >= 0
            except ValueError:
                print "All inputs must be real numbers greater than or equal to zero"
                
        return q*s*cd
    
    def lift(self,q,s,cl):
        
        inputs = [q,s,cl]
        
        
        for x in inputs:
            try:
                x = float(x)
                x >= 0
            except ValueError:
                print "All inputs must be real numbers greater than or equal to zero"
                
        return q*s*cl
    
            
            
            
