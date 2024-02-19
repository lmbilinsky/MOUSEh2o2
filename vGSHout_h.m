
%  GSH export function
% high affinity sinusoidal export

function a = vGSHout_h(gsh_cyt);

km = 150; 
vmax=306;

a = vmax * gsh_cyt./(km + gsh_cyt);

