import pandas as pd
from moi_calc_calculator import *

'''
NOTE:
    always assumes that that orientation of cylinders and hoops 
    face the Z axis; i.e. the axis of symmetry is always assumed to be 
    the z axis.
'''

#reading in data
#This filename needs to be changed when working with other mass budgets,
#and depending on folder opened.
file_name = 'MOI_Calc/Ouroboros Mass Budget - Ouroboros Mass Budget (1).csv'

rawData = pd.read_csv(file_name, sep=',', dtype=str)

#getting header
header = rawData.head(0)


def getPartsNeedingMOICalc() -> dict:
    '''
    Returns dict of all parts with a valid shape column from
    the Mass budget as key with value being the matching calculation.
    '''
    #{shape : calculation}
    shapesAndCalcs = {
        "Cylindrical" : CalculateMOI.calculate_cylinder,
        "Hoop" : CalculateMOI.calculate_hoop,
        "Spherical" : CalculateMOI.calculate_spherical
    }
    
    partsWithCalculation = {}
    indicesWithCalculation = {}
    
    #Gather useful data and fill the dictionary
    #The dictionary will contain useful data's indicies
    #as the key and the matching calculation function
    #as the value
    i = 0
    while i < rawData['class'].size:
        if (rawData['Shape'][i] != 'nan') and (type(rawData['Shape'][i]) != float):
            indicesWithCalculation[i] = shapesAndCalcs[rawData['Shape'][i]]
            partsWithCalculation[rawData['Part'][i]] = shapesAndCalcs[rawData['Shape'][i]]
        i += 1
    
    partsIndicesCalc = [partsWithCalculation, indicesWithCalculation]
    return indicesWithCalculation

def calculateMomentsOfInertia() -> dict:
    '''
    Returns a dictionary with a part name as key
    and list of moments of inertia in X, Y, and Z as value
    '''
    partindices = getPartsNeedingMOICalc()
    indices = [key for key in partindices]
    
    dictOfMOIs = {}
    
    for index in indices:
        #getting parameters
        vector = []

        if partindices[index] == CalculateMOI.calculate_cylinder:
            mass = float(rawData['kg'][index])
            length = float(rawData['length'][index])
            
            if (rawData['outerDiameter'][index] != 'nan') and (type(rawData['outerDiameter'][index]) != float):
                outer_diameter = float(rawData['outerDiameter'][index])
                #calculate here
                vector.append(partindices[index](mass, outer_diameter, length, axis = 'x'))
                vector.append(partindices[index](mass, outer_diameter, length, axis = 'y'))
                vector.append(partindices[index](mass, outer_diameter, length, axis = 'z'))
                
            elif (rawData['startDiameter'][index] != 'nan') and (type(rawData['startDiameter'][index]) != float):
                if (rawData['endDiameter'][index] != 'nan') and (type(rawData['endDiameter'][index]) != float):
                    outer_diameter = float(rawData['startDiameter'][index])
                    inner_diameter = float(rawData['endDiameter'][index])
                    #calculate
                    vector.append(partindices[index](mass, outer_diameter, length, inner_radius=inner_diameter, axis = "x"))
                    vector.append(partindices[index](mass, outer_diameter, length, inner_radius=inner_diameter, axis = "y"))
                    vector.append(partindices[index](mass, outer_diameter, length, inner_radius=inner_diameter, axis = "z"))
                    
                    
        elif partindices[index] == CalculateMOI.calculate_hoop:
            mass = float(rawData['kg'][index])
            
            if (rawData['outerDiameter'][index] != 'nan') and (type(rawData['outerDiameter'][index]) != float):
                outer_diameter = float(rawData['outerDiameter'][index])
                #calculate
                vector.append(partindices[index](mass, outer_diameter, axis = 'x'))
                vector.append(partindices[index](mass, outer_diameter, axis = 'y'))
                vector.append(partindices[index](mass, outer_diameter, axis = 'z'))
                
        elif partindices[index] == CalculateMOI.calculate_spherical:
            mass = float(rawData['kg'][index])
            
            if (rawData['outerDiameter'][index] != 'nan') and (type(rawData['outerDiameter'][index]) != float):
                outer_diameter = float(rawData['outerDiameter'][index])
                #calcualte
                vector.append(partindices[index](mass, outer_diameter, type = 'solid'))
        
        dictOfMOIs[rawData['Part'][index]] = vector
            
    return dictOfMOIs

#Set this to add directly to the mass budget. 
finalCalcs = calculateMomentsOfInertia()

for n in finalCalcs:
    print(f'{n} : {finalCalcs[n]}')

    
#First Sucessful test run  
#These are nonsense results as the shapes don't always match
#The given 
'''
Nose Cone : [0.060639479077532, 0.060639479077532, 0.019697625]
Nose Cone Tip Bolt : []
Nose Cone Tip : [8.602273376064e-05]
Nose Cone Shell : [0.05321744187795887, 0.05321744187795887, 0.014365134663221597]
Chutes : [nan, nan, 0.010118875000000001]
Recovery Tube : [0.036515820765665855, 0.036515820765665855, 0.011226375]
Avionics Tube : [0.036515820765665855, 0.036515820765665855, 0.011226375]
Payload Interference Plate : []
Payload Inserts : []
Pressure Vessel : [8.051331931673106, 8.051331931673106, 0.29005263939033515]
Upper End Cap : [0.0032239608108651603, 0.0032239608108651603, 0.0063661852514750004]
Lower End Cap : [0.0032239608108651603, 0.0032239608108651603, 0.0063661852514750004]
Upper Casing : [0.9335222061983164, 0.9335222061983164, 0.9476144125480296]
Lower Casing : [0.8623019258521967, 0.8623019258521967, 0.027516614046210225]
Bolts : []
Tube : [3.0540599947075133, 3.0540599947075133, 0.12083737499999998]
Upper O-rings : [nan, nan, 4.212538687809e-05]
Lower O-rings : [nan, nan, 4.212538687809e-05]
Upper : []
Lower : []
Fins : []
Fin Can : [0.17799629579329754, 0.17799629579329754, 0.039693294226672136]
Fuel Grain casing With Mixing Plate : [nan, nan, nan]
Chamber : [0.1184602026159, 0.1184602026159, 0.0374085]
Nozzle : [0.009972910503675004, 0.009972910503675004, 0.017241883838744004]
Nozzle Housing : [nan, nan, 0.05955067964007994]
Bolt Ring 1 : [nan, nan, 0.004445280000000001]
Injector Bulkhead : [0.007237292105081999, 0.007237292105081999, 0.01378275]
Bolt Ring 2 : [nan, nan, 0.006975475000000001]
'''
