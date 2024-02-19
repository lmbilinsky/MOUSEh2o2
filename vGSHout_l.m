% GSH export function
% low affinity sinusoidal export

function a = vGSHout_l(gsh_cyt);

km = 3000; 

%Old info.: efflux from liver of rat is given in Fig. 2, Ookhtens85: 15 nm/(min*gram). Bartoli78 states a slightly smaller efflux. We don't use this info. in the mouse model. The vmax we choose for mouse is to be in agreement with the half-life of GSH in mouse. Mouse has a higher liver concentration of GSH than rat, so efflux must be larger.

vmax=2247;

a = vmax * (gsh_cyt.*gsh_cyt.*gsh_cyt)./ (km^3 + gsh_cyt.*gsh_cyt.*gsh_cyt);

