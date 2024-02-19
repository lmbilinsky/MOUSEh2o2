Update to the model: for increased realism, the per-cell death rate of hepatocytes is now 0 if the intracellular h2o2 concentration, [h2o2], is less than 1 micromolar. Otherwise, it is as described in the paper

To simulate the in vivo computational model of endogenous H2O2 metabolism in mouse liver, run MOUSEh2o2.m. All other .m files (which calculate reaction and transport velocities) in this repository are called by 
this master function and should be stored in the same directory/folder as it. The file FinalState.txt lists the final concentrations and final reaction and transport velocities at the end of the simulation
period.
