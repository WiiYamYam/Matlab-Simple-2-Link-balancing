# Matlab-Simple-2-Link-balancing
This project was for learning the control systems required to balance a simple 2 link legged balancing robot. It is not to design a robot for manufacturing but to practise the math that makes the robot balance. This simulation does not deal with the non-linear model as this is beyond my understanding right now. Perhaps in later projects things like this will be present. The purpose of this project was also to document my first successful steps in dynamics of robotic systems.


included are 6 files

vvvvvvvvvvvvvvv
model.slx - Matlab simscape multibody simulation of the system. You can play with the block diagram here.

generate_lqr_data.m - This file generates LQR K gain values for a series of L0's and saves them to a lookup table for later use.

get_K_L0.m - This file read the robot's current L0 (length of the imaginary line between the hip and ankle) and finds the LQR K gains that are a closest match.

leg_conv.m - This file does the forward kinematics of the simple 2 link planar chain that is the leg. It then uses jacobians to calculate the required torque to drive the knee into the correct position to change the hight of L0. It then uses symbolic toolbox to export this to a function that doesnt have any functions that simulink doesnt like.
leg_convolution.m - simulink doesnt liek using some matlab functions so this is the function called to solve for the torques in the VMC.

K_lookup_data_good.mat - this is the data saved by generate_lqr_data.m. In it are various K gain matrices and their corresponding L0.
