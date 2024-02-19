%catalase removal of H2O2

function a = Vcat(h2o2);




%Vmax=1e6;    %science direct: Heme Proteins, J. Everse, Encyclopedia of Biological Chemistry 2013 

Vmax=1e6;



%Everse says that under ideal conditions each subunit can dismute 2e5 mol
%H2O2 per second

km = 32300;   %for mouse liver 


a = Vmax*h2o2./(km + h2o2);  %Chen 2015, for mouse liver



