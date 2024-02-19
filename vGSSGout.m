%  GSSG export

function a = vGSSGout(gssg_cyt);

km = 400;  


vmax=180;


a = vmax * gssg_cyt./(km + gssg_cyt);

