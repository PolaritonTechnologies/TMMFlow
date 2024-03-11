# Author: Florian Le Roux
import numpy as np

class CalculationInfo():
    
#this class prepares all the info necessary for calculation
    def __init__(self, intransmission, angleMin, angleMax, angleStep, wavelengthMin, wavelengthMax, wavelengthStep, structure, polarization, azimuthalAngles, inEnergy, Zmin, Zmax):
    
    #can add incident medium and substrate, for now default to air and glass
        self.intransmission = intransmission
        self.angleMin = angleMin
        self.angleMax = angleMax
        self.angleStep = angleStep
        self.wavelengthMin = wavelengthMin
        self.wavelengthMax = wavelengthMax
        self.wavelengthStep = wavelengthStep
        self.structure = structure   
        self.polarization = polarization
        self.azimuthalAngles = azimuthalAngles
        self.calculationOrderName = self.formCalculationName()
        self.inEnergy = inEnergy
        self.EMin = 1239.82/float(wavelengthMax)
        self.EMax = 1239.82/float(wavelengthMin)
        self.Zmin = Zmin
        self.Zmax = Zmax

    def formCalculationName(self):
    
        structureAsString = ""
        
        if self.azimuthalAngles is np.array:  
            
            azimuthalMin = self.azimuthalAngles[0]
            azimuthalMax = self.azimuthalAngles[np.size(self.azimuthalAngles)-1]         
            
            for element in self.structure:
                
                structureAsString = structureAsString + str(element[0]) + str(element[1])     
            
            return str(structureAsString) + str(self.angleMin) + str(self.angleMax) + str(self.angleStep) + str(self.wavelengthMin) + str(self.wavelengthMax) + str(self.wavelengthStep) + 'A' + str(azimuthalMin) + str(azimuthalMax) + '.txt'            
        
        else:
            
            aziMin = self.azimuthalAngles
            
            for element in self.structure:
                
                structureAsString = structureAsString + str(element[0]) + str(element[1])     
            
            return str(structureAsString) + str(self.angleMin) + str(self.angleMax) + str(self.angleStep) + str(self.wavelengthMin) + str(self.wavelengthMax) + str(self.wavelengthStep) + 'AtA' + str(aziMin) + '.txt'            
                    