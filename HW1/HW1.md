The objective of this homework set is to help you understand the first steps to simulating dynamics systems.

1.  Carefully dissect hw02_armSim.py, and understand what every line of code in that file is doing.

For each homework problem, you will be given:
a parameter file like armParam.py that contains the physical parameters
an animation file like armAnimiation.py that draws a graphic
 dataPlotter.py that plots the state and input data
Study these files and try to understand how they work.

- [x] Done, I feel like I understand how they work.

2.  In a similar way, carefully dissect hw02_pendulumSim.py and hw02_satelliteSim.py.  Run these launch files and understand what you are observing.

 

The hw02_*.py files simply plot the systems but do not implement the dynamics.  The dynamics will be implemented in the next homework set.

- [x] Done, I feel like I understand how the code works, not the EOMs though... 

3.  In the _D_mass directory, create the simulation launch file hw02_massSim.py, and use it to animate the mass moving back and forth according to a sin wave with amplitude of 1 (m), frequency of 0.5 Hz, and y_offset of 0.2 (m).

 

3.  In the _E_blockBeam directory, create the simulation launch file hw02_blockBeamSim.py, and use it to animate the block-on-beam system where the block is moving along the beam according to a sin wave with amplitude of 0.05 (m), frequency of 0.5 Hz, and y_offset of 0.2 (m), and where the angle of the beam is following a square wave with amplitude pi/8 (rad), and frequency of 0.1 Hz.

 

4.  In the _F_VTOL directory, create the simulation launch file hw02_VTOLSim.py, and use it to animate the planar VTOL system where the the position z follows a sin wave with amplitude of 4 (m), frequency of 0.1 Hz, and y_offset of 5 (m), altitude h follows a square wave with amplitude of 2 (m), frequency of 0.1 Hz, and