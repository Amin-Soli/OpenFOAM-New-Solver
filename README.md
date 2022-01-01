
# This solver is called "simpleFoamWithPhaseChangeTwoPhaseMixture", and only executable in OpenFOAM 4.0.

# This solver solves two phase mixture equations in which phase change occurs. 

# This solver is steady state and uses two different models called Schnerr-Sauer and Zwart in order to calculate mass exchange rate.

# In order to have better stability, the volume fraction (alpha) equation is solved for vapour instead of water. Besides, a diffusion term is added to alpha equation that makes the solver have better stability.

# You can control the diffusion term by a coefficient defined in "transportProperties" file in "constant" folder as "alfaCoeffDiffusion".

# This solver is also able to calculate forces and torques on patches for two phase mixture.

# There is also a "commented section" at the end of "simpleFoamWithPhaseChangeTwoPhaseMixture.C" file in which the solver calculates forces and torques on a specific patch for two phase mixture. If you want to use this feature, please uncomment that.

# If you want to compile this solver, all you need to do is:
  1. Open terminal
  2. Go to "solver/Make" address
  3. First type "wclean" and then type "wmake" in the terminal

# Note: In this solver, water represents vapour ,and vapour represents water.

# There is also a case study in which you can test the solver.

# In order to run this case study, all you need to do is:
  1. Open terminal
  2. go to "caseStudy" address
  3. type "./Allrun" in the terminal

# Note: alpha.water represents vapour volume fraction.

