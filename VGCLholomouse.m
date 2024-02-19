% Glutamate cysteine ligase (GCL) 
% from Chen05

function a = VGCLholomouse(cys,gsh_cyt);

glu=1400;
atp=3000;


V=235820;  %chosen to give agreement with Sekura74's value of 145 minutes for the half-life of GSH in mouse liver  (title of paper: "GSH turnover in kidney")



k1 = 220;       %Km cys
k2 = 480;       %Km glut
k3= 870;        %km ATP
K2is=800;
K2ii=3100;
K3is=6500;
K3ii=3900;

    

a = V*(cys./(k1+cys)).*(glu./(k2+k2.*gsh_cyt./K2is+glu+glu.*gsh_cyt./K2ii)).*(atp./(k3+k3.*gsh_cyt./K3is+atp+atp.*gsh_cyt./K3ii));  

%a = V*(cys./(k1+cys)).*(glu./(k2+k2.*gsh_cyt./K2is+glu+glu.*gsh_cyt./K2ii));

%.*(atp./(k3+k3.*gsh_cyt./K3is+atp+atp.*gsh_cyt./K3ii));  

