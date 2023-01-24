The objective of this file is to implement the dynamics, or physics engine for the mass, the blockBeam, and the VTOL.

 

1.   Carefully dissect hw03_armSim.py, making sure that you understand every line of code, especially how the dynamics are called in the systems.  Similarly study  hw03_pendulumSim.py, and hw03_satelliteSim.py.

  - [ ]  

2.  Carefully study armDynamics.py, again making sure that you understand every line of code.  In this case you may want to look at the dynamic equations of motion for the arm, to see how the differential equation is represented and solved.  Do the same thing for pendulumDynamics.py and satelliteDynamics.py.

 

3.  In the _D_mass directory, create the file massDynamics.py and implement the dynamics equations of motion for the mass-spring-damper system.  Use the file test_dynamics.py to test your dynamics file to ensure that it is providing the same outputs as my files for different initial conditions and inputs.  Create the file hw03_massSim.py, that simulates the mass system when the input force is a sin wave with amplitude 10 (N), and frequency of 1 Hz.  Note:  The simulation will not do anything interesting.

 - [ ] ***mostly base this off of the arm file, the vtol base off of pendulum or satellite, the equations are in the slides! Take in a state and a u, the spit out an xdot, the derivative of each of the states

4.  In the _E_blockBeam directory, create the file blockBeamDynamics.py and implement the dynamics equations of motion for the blockBeam system.  Use the file test_dynamics.py to test your dynamics file to ensure that it is providing the same outputs as my files for different initial conditions and inputs.  Create the file hw03_blockBeamSim.py, that simulates the blockBeam system when the input force is a sin wave with amplitude 0.5 (N), and frequency of 1 Hz, and a y_offset of 11.5 (N).  Note:  The simulation will be unstable.

 

5.  In the _F_VTOL directory, create the file VTOLDynamics.py and implement the dynamics equations of motion for the planar VTOL system.  Use the file test_dynamics.py to test your dynamics file to ensure that it is providing the same outputs as my files for different initial conditions and inputs.  Create the file hw03_VTOLSim.py, that simulates the planar VTOL system when the thrust on both motors is a sin wave with amplitude 10 (N), and frequency of 1 Hz.  Note:  The simulation will not do anything interesting.