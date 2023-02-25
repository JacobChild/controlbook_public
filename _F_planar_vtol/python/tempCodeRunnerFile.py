#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlPIDhw10 import ctrlPD

#instantiate mass, controller, and reference classes
VTOL = VTOLDynamics(alpha = 0.2) #Instantiated the blockbeamDynamics class as blockbeam in this file
controller = ctrlPD()
#instantiate reference input classes
zPosRefSig = signalGenerator(3.0, 0.08, 2.5) #amplitude, frequency, y_offset
hInputRefSig = signalGenerator(5.0, 0.01, 5.0) #this will be where the torque input/controller is definded

#instantiate the simulation plots and animation