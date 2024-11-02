import pandas as pd
import math 

'''
NOTE:
    always assumes that that orientation of cylinders and hoops 
    face the Z axis; i.e. the axis of symmetry is always assumed to be 
    the z axis.
'''

#TypeError class for spherical MOI calculation
class TypeError(Exception):
    def __init__(self, message = "argument type should be solid or thin"):
        self.message = message
        super().__init__(self.message)
        
#AxisError for cylinders and hoops 
class AxisError(Exception):
    def __init__(self, message = "argument axis should be z, x or y"):
        self.message = message
        super().__init__(self.message)

class CalculateMOI:
    def calculate_cylinder(mass:float, outer_radius:float, length:float, inner_radius:float = 0, axis:str = "z") -> float:
        '''
        takes in attributes of a cylindical or disc component and returns the moment of inertia about
        a chosen axis
        
        Args:
            mass of component
            outer radius of component - MUST be >0
            inner radius of component - should be 0 if the cylinder is solid
            length of component
            axis of the moment: either z, x or y
        returns:
            moment of inertia about chosen axis (default z)
            
        all returns from known formulas 
        '''
        #pre-check that numbers are valid
        if mass < 0 or outer_radius < 0 or length < 0 or inner_radius < 0:
            raise ValueError
        
        #calculate about z axis
        if axis == "z":
            return 0.5*(mass)*(outer_radius**2 + inner_radius**2)
        
        #calculate about xy
        elif axis == "x" or axis == "y":
            return (0.0833333)*(mass)*(3*(inner_radius**2 + outer_radius**2) + (length**2))
        
        #incase bad inputs
        else:
            print("axis should be symmetry or diamter")
            raise TypeError
        
            
    def calculate_hoop(mass:float, radius:float, axis:str = "symmetry") -> float:
        '''
        used for infinitely thin hoops 
        axis can be z, x or y.
        '''
        
        #about z axis
        if axis == "z":
            return mass*(radius**2)
        
        #about xy
        elif axis == "x" or axis == "y":
            return 0.5*(mass)*(radius**2)
        
        #error check
        else:
            print("axis should be symmetry or diamter")
            raise TypeError
        
        
    def calculate_spherical(mass:float, radius:float, type:str = "solid") -> float:
        '''
        similar args to previous functions, with type of sphere instead of axis 
        type should be solid or thin 
        '''
        
        #if solid
        if type == "solid":
            return (0.4)*(radius**2)*(mass)
        
        #if thin
        elif type == "thin":
            return (0.4)*(radius**2)*(mass)
        
        #bad input
        else:
            print("type of sphere should be solid or thin")
            raise TypeError
        
    



