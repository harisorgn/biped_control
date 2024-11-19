# Bipedal Walking Control

## Introduction
This project includes the code that was developed to simulate leg placement strategies in a 3D spring-mass bipedal model. These strategies rely on the body dynamics of the 3D bipedal model and do not contain any feedforward or feedback components. Thus, they are energetically efficient and easier to implement in a physical bipedal robot. Furthermore, they are inspired by how humans rely on the elasticity of their joints and tendons while walking, so that they do not have to adjust for any minor perturbation to their walking trajectory. 


The package implements the bipedal model as a hybrid dynamical system that switches between a single-leg and a double-leg phase, a custom root-finding method to identify limit cycles of steady walking, utility functions to filter the found walking configurations given user's criteria, measurement functions that calculate the ground forces that apply on each leg and plotting functions to visualize different parts of the analysis.

There are currently two successful strategies that have been identified, `ratio-beta` and `beta`. The analysis code for each one is in the respective folder.

This work was a dissertation project for the MSc Robotics programme (2016-2017) at the University of Bristol and the report titled “Leg placement control based on passive dynamics of a 3D walking model” is included as well for a detailed description of the applied methods.


## How to use
The functions in both folders (leg strategies) are identical for the most part, with some strategy-specific parts noted with comments. Different strategies can be applied to the model and their stability assessed by changing these blocks of code in `biped3_step` and `biped3_TD` functions, which include the step dynamics and the touch-down event respectively. Some additions may be necessary in the systems structure sys if novel control variables are introduced.

The sequence of the code is as follows :
- First the `biped3_lc_search` function, which utilises a custom-built Newton-Raphson algorithm to find stable limit cycles of the system. 
- Then `biped3_lc_param` parses through the found limit cycles, identifies any potential duplicate walking patterns that correspond to these limit cycles and calculates useful gait parameters, such as step time, length, width and mean velocities. 
- From that point on the `*_lc.mat` file (where * is either `ratio_beta` or `beta`, depending on the considered leg strategy) containing all the kinematic and control parameters of each unique limit cycle can be used for further analysis. 

## Examples
Some examples of using the current package are presented in functions 
- `biped3_step_ex` to simulate walking
- `biped3_grf_ex` to calculate and plot ground reaction forces 
- `biped3_basin_ex` to calculate partial basins of attraction for a limit cycle 
- `biped3_basin_all` to calculates partial basins of attraction for all found limit cycles. 
- `biped3_basin_plot` is used to plot the resulted basins of attraction, either for a single pattern or for all of them.

Developed by Haris Organtzidis
