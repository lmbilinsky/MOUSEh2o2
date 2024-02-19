function MOUSEh2o2


load yinitialcondition

%GSH MODEL PARAMETERS
Vplasmatotal=0.00113943;

NADPH = 50;
diffh2o2=2284581;
vh2o2other=3026;
vh2o2prod=28567;
gsh_deg=17.62;

cellnum=1.8e8;

%%%%%%%%%%%%% FUNDAMENTAL CONSTANTS
hematocrit=0.45;
BW=0.025; 
QCC=16.5;


VbloodC=0.049;
VliverC=0.055;
VkidneyC=0.017;
VlungC=0.007;         %
VurinarybladderC=0.0009;
VslowlyC=0.82;        %sum of fat and slowly perfused in Andersen et al. 1987
VrichlyC=0.09;        %sum of liver and richly perfused in Andersen et al. 1987

VcapliverC=0.31;
VtissliverC=0.69;
VcapkidneyC=0.24;
VtisskidneyC=0.76;
VcaplungC=0.5;
VtisslungC=0.5;
VcapurinarybladderC=0.03;
VtissurinarybladderC=0.97;
VcapslowlyC=0.03;  %Value for skin from Brown97
VtissslowlyC=0.97;
VcaprichlyC=0.31; %Value for liver from Brown97
VtissrichlyC=0.69;

QliverC=0.161;
QkidneyC=0.091;
QlungC=0.005;         
QurinarybladderC=0.0033;
QslowlyC=0.24;        %sum of fat and slowly perfused in Andersen et al. 1987  
QrichlyC=0.76;        %calculated 



%%%%%%%%%%%%% COMPUTED CONSTANTS
QCP=(1-hematocrit)*QCC*BW^0.75;
Vblood=VbloodC*BW;

Vliver=VliverC*BW;
Vcapliver=VcapliverC*Vliver;
Vtissliver=VtissliverC*Vliver;
Vkidney=VkidneyC*BW;
Vcapkidney=VcapkidneyC*Vkidney;
Vtisskidney=VtisskidneyC*Vkidney;
Vlung=VlungC*BW;
Vcaplung=VcaplungC*Vlung;
Vtisslung=VtisslungC*Vlung;
Vurinarybladder=VurinarybladderC*BW;
Vcapurinarybladder=VcapurinarybladderC*Vurinarybladder;
Vtissurinarybladder=VtissurinarybladderC*Vurinarybladder;
Vrichly=VrichlyC*BW - Vliver - Vkidney - Vlung - Vurinarybladder;
Vcaprichly=VcaprichlyC*Vrichly;
Vtissrichly=VtissrichlyC*Vrichly;
Vslowly=VslowlyC*BW - Vblood;
Vcapslowly=VcapslowlyC*Vslowly;
Vtissslowly=VtissslowlyC*Vslowly;

Qliver=QliverC*QCP;
Qkidney=QkidneyC*QCP;
Qlung=QlungC*QCP;         
Qurinarybladder=QurinarybladderC*QCP;
Qslowly=QslowlyC*QCP;        %from Fisher00 (trichloroethylene).
Qrichly=QrichlyC*QCP - Qliver -Qkidney - Qlung - Qurinarybladder;        %from Fisher00 







yinit_noxenobiotic=zeros(265,1);
yinit_noxenobiotic(201:265)=yinit;
[t_highdose y_highdose]=ode15s(@mouseh2o2,[0 24],yinit_noxenobiotic);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE SHOWING WHAT'S GOING ON INSIDE THE LIVER
living=y_highdose(:,265)*100/cellnum;
h2o2=0.85*y_highdose(:,232)+0.15*y_highdose(:,221);
gsh=0.85*y_highdose(:,205)+0.15*y_highdose(:,203);
gssg=0.85*y_highdose(:,212)+0.15*y_highdose(:,230);


figure 
subplot(4,1,1)
plot(t_highdose,h2o2,'Linewidth',2)
xlim([0 24])
title('H2O2 IN LIVER CELLS (microM)')
subplot(4,1,2)
plot(t_highdose,gsh,'Linewidth',2)
xlim([0 24])
title('GSH IN LIVER')
subplot(4,1,3)
plot(t_highdose,gssg,'Linewidth',2)
xlim([0 24])
title('GSSG IN LIVER')
subplot(4,1,4)
plot(t_highdose,y_highdose(:,265)*100/cellnum,'Linewidth',2)
xlim([0 24])
title('LIVING HEPATOCYTES')











%%%%%%%%%%%%%%%%%%%%%%% PRINT-OUT OF FINAL STATE
T=t_highdose;
Y=y_highdose;

Vcapliverplasma=Vcapliver*(1-hematocrit);

%AA INPUT INTO PLASMA
cysinput=cysin(T(length(T)));


%CONCENTRATIONS IN PLASMA
cystine_plasma=Y(length(T),241);
cys_plasma = Y(length(T),242);
gsh_plasma = Y(length(T),243);


%CONCENTRATIONS IN LIVER CAPILLARY PLASMA
cys_lcp = Y(length(T),244);
gsh_lcp = Y(length(T),245);


%CONCENTRATIONS IN LIVER
cys_mz  = Y(length(T),204);
glutamylcys_mz = Y(length(T),206);
gsh_cyt_mz  = Y(length(T),205);
gsh_mito_mz=Y(length(T),203);
gssg_cyt_mz = Y(length(T),212);
gssg_mito_mz=Y(length(T),230);
h2o2_cyt_mz=Y(length(T),232);
h2o2_mito_mz=Y(length(T),221);


%VELOCITIES
vlcpCYSl=Vcysin(cys_lcp);
vlGSHlcp=vGSHout_l(gsh_cyt_mz) + vGSHout_h(gsh_cyt_mz); 
vplasmaGSHdeg=gsh_deg*gsh_plasma;
vlcpCYSplasma=cys_lcp*Qliver/Vcapliverplasma;
vlcpGSHplasma=gsh_lcp*Qliver/Vcapliverplasma;
vplasmaCYSlcp=cys_plasma*Qliver/Vplasmatotal;
vplasmaGSHlcp=gsh_plasma*Qliver/Vplasmatotal;


vcysdioxy_mz=Vcysdioxygenase(cys_mz);
vgcl_mz=  VGCLholomouse(cys_mz,gsh_cyt_mz);
vgss_mz = VGSSrat(glutamylcys_mz);
vgshcyt_to_mito_mz=Vgshcyttomito(gsh_cyt_mz);
vgshmito_to_cyt_mz=Vgshmitotocyt(gsh_mito_mz);
vdiffh2o2mito_to_cyt_mz=diffh2o2*h2o2_mito_mz;
vcat_mzcyt=Vcat(h2o2_cyt_mz);
vgpx_mzcyt=VGPXcyt(gsh_cyt_mz,h2o2_cyt_mz);
vgr_mzcyt=VGRcyt(gssg_cyt_mz,NADPH);
vgpx_mzmito=VGPXmito(gsh_mito_mz,h2o2_mito_mz);
vgr_mzmito=VGRmito(gssg_mito_mz,NADPH);
vmzGSSGbile= vGSSGout(gssg_cyt_mz);

M=zeros(116,1);
M(4)=cys_mz;
M(5) = gsh_cyt_mz;
M(6) = glutamylcys_mz;
M(8) = 0;
M(10)=0;
M(11)=0;
M(12)=gssg_cyt_mz;
M(17)=vcysdioxy_mz;
M(19)=vgcl_mz;
M(20)=vgss_mz;
M(27)=vlcpCYSl;
M(28)=vlGSHlcp;
M(29)=vmzGSSGbile;
M(30)=vgpx_mzcyt;
M(31)=vgr_mzcyt;
M(33)=vplasmaGSHdeg;
M(51)=vgshcyt_to_mito_mz;
M(53)=vgpx_mzmito;
M(54)=vgr_mzmito;
M(55)=vgshmito_to_cyt_mz;
M(57)=gsh_mito_mz;
M(58)=gssg_mito_mz;
M(89)=vh2o2prod;
M(91)=h2o2_mito_mz;
M(98)=vdiffh2o2mito_to_cyt_mz;
M(99)=h2o2_cyt_mz;
M(103)=vh2o2other;
M(104)=vlcpCYSplasma;
M(105)=vlcpGSHplasma;
M(106)=vplasmaCYSlcp;
M(107)=vplasmaGSHlcp;
M(108)=cystine_plasma;
M(109)=cys_plasma;
M(110)=gsh_plasma;
M(111)=cys_lcp;
M(112)=gsh_lcp;
M(113)=vcat_mzcyt;



PrintFinalValuesMousePBPK(Vplasmatotal, Vtissliver, Vcapliverplasma, M, cysinput);


end




function dydt=mouseh2o2(t,y)

%GSH MODEL PARAMETERS
Vplasmatotal=0.00113943;

NADPH = 50;
diffh2o2=2284581;
vh2o2other=3026;
vh2o2prod=28567;
gsh_deg=17.62;

cellnum=1.8e8;

death=0.0009;
h2o2=0.85*y(232)+0.15*y(221);


limit=150;
if h2o2 < 1
     deathrate=0;
elseif ((h2o2 >= 1) & (h2o2 < limit))
     deathrate=death*h2o2;
else
     deathrate=death*limit;
end



%limit=150;
%if h2o2<limit
%    deathrate=death*h2o2;
%else
%    deathrate=death*limit;
%end




%%%%%%%%%%%%% FUNDAMENTAL CONSTANTS
hematocrit=0.45;
BW=0.025; 
QCC=16.5;


VbloodC=0.049;
VliverC=0.055;
VkidneyC=0.017;
VlungC=0.007;         %
VurinarybladderC=0.0009;
VslowlyC=0.82;        %sum of fat and slowly perfused in Andersen et al. 1987
VrichlyC=0.09;        %sum of liver and richly perfused in Andersen et al. 1987

VcapliverC=0.31;
VtissliverC=0.69;
VcapkidneyC=0.24;
VtisskidneyC=0.76;
VcaplungC=0.5;
VtisslungC=0.5;
VcapurinarybladderC=0.03;
VtissurinarybladderC=0.97;
VcapslowlyC=0.03;  %Value for skin from Brown97
VtissslowlyC=0.97;
VcaprichlyC=0.31; %Value for liver from Brown97
VtissrichlyC=0.69;

QliverC=0.161;
QkidneyC=0.091;
QlungC=0.005;         
QurinarybladderC=0.0033;
QslowlyC=0.24;        %sum of fat and slowly perfused in Andersen et al. 1987  
QrichlyC=0.76;        %calculated 

PIdma3slowlyC=3;
PIdma3richlyC=1;
PIdma3lungC=0.5;
PIdma3kidneyC=2;
PIdma3liverC=0.25;
PIdma3urinarybladderC=0.01;

PIdma5slowlyC=0.3;
PIdma5richlyC=0.1;
PIdma5lungC=0.001;
PIdma5kidneyC=0.2;
PIdma5liverC=0.02;
PIdma5urinarybladderC=0.001;

Pdma3slowly=500;
Pdma3richly=200;
Pdma3lung=200;
Pdma3kidney=600;
Pdma3liver=100;
Pdma3urinarybladder=300;

Pdma5slowly=5;
Pdma5richly=2;
Pdma5lung=2;
Pdma5kidney=6;
Pdma5liver=1;
Pdma5urinarybladder=3;

dma3todma5liverC=0.1;
dma3todma5kidneyC=0.03;
dma3todma5urinarybladderC=0.002;
dma3todma5lungC=0.01;
dma3todma5slowlyC=2;
dma3todma5richlyC=0.02;
dma3todma5rbcC=0.04;

plasmatorbc5C=0.008;
rbctoplasma5C=0.02;
plasmatorbc3C=0.2;
rbctoplasma3C=0.03;

dma5absorpKmC=2000;
dma5reductiongutKmC=800;
dma3absorpC=1;
dma5absorpVmaxC=50;
dma5reductiongutVmaxC=120;
dma3totmaocolonC=1;
tocolonC=0.25;
dmatofecesC=0.4;
dma3tourineC=1;
dma5tourineC=1;
tmaotourineC=0.3;


%%%%%%%%%%%%% COMPUTED CONSTANTS

QCP=(1-hematocrit)*QCC*BW^0.75;
Vblood=VbloodC*BW;
Vliver=VliverC*BW;
Vcapliver=VcapliverC*Vliver;
Vtissliver=VtissliverC*Vliver;
Vkidney=VkidneyC*BW;
Vcapkidney=VcapkidneyC*Vkidney;
Vtisskidney=VtisskidneyC*Vkidney;
Vlung=VlungC*BW;
Vcaplung=VcaplungC*Vlung;
Vtisslung=VtisslungC*Vlung;
Vurinarybladder=VurinarybladderC*BW;
Vcapurinarybladder=VcapurinarybladderC*Vurinarybladder;
Vtissurinarybladder=VtissurinarybladderC*Vurinarybladder;
Vrichly=VrichlyC*BW - Vliver - Vkidney - Vlung - Vurinarybladder;
Vcaprichly=VcaprichlyC*Vrichly;
Vtissrichly=VtissrichlyC*Vrichly;
Vslowly=VslowlyC*BW - Vblood;
Vcapslowly=VcapslowlyC*Vslowly;
Vtissslowly=VtissslowlyC*Vslowly;



Qliver=QliverC*QCP;
Qkidney=QkidneyC*QCP;
Qlung=QlungC*QCP;         
Qurinarybladder=QurinarybladderC*QCP;
Qslowly=QslowlyC*QCP;        %from Fisher00 (trichloroethylene).
Qrichly=QrichlyC*QCP - Qliver -Qkidney - Qlung - Qurinarybladder;        %from Fisher00 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE STATE VARIABLE FLUX
fracliver=y(265)/cellnum;

dydt=zeros(265,1);


dydt(203)=Vgshcyttomito(y(205))*(0.85/0.15)-Vgshmitotocyt(y(203)) - 2*VGPXmito(y(203), y(221)) + 2*VGRmito(y(230), NADPH); %liver mito GSH

dydt(230)= VGPXmito(y(203),y(221)) - VGRmito(y(230), NADPH); % liver mito GSSG 

dydt(204)=((Vcapliver*(1-hematocrit))/Vtissliver)*Vcysin(y(244)) - VGCLholomouse(y(204),y(205)) + 100 - Vcysdioxygenase(y(204)); % liver cysteine. Extra 100 is input from MET cycle

dydt(205)=VGSSrat(y(206)) - 2*VGPXcyt(y(205), y(232)) + 2*VGRcyt(y(212), NADPH) - vGSHout_h(y(205)) - vGSHout_l(y(205)) - Vgshcyttomito(y(205)) + Vgshmitotocyt(y(203))*(0.15/0.85); %liver cytosolic GSH. 

dydt(206)=VGCLholomouse(y(204),y(205))/0.85 - VGSSrat(y(206)); %liver gamma-GC. 

dydt(212)=VGPXcyt(y(205),y(232)) - VGRcyt(y(212), NADPH) - vGSSGout(y(212)); % liver cytosolic GSSG

dydt(221)= vh2o2prod - diffh2o2*y(221) - VGPXmito(y(203),y(221)); %H2O2 in liver mitochondria. 

dydt(232)= (0.15/0.85)*diffh2o2*y(221) - VGPXcyt(y(205), y(232)) + vh2o2other  - Vcat(y(232)); %H2O2 in liver cytosol




dydt(257) = sin(100*t);

dydt(265) = - deathrate*y(265); % living hepatocytes (neglect regeneration)

dydt(241) = 100*y(242) - 100*(24/50)*y(241); %plasma cystine

dydt(242)= cysin(t) + gsh_deg*y(243) - (Qliver/Vplasmatotal)*y(242) + (Qliver/Vplasmatotal)*y(244) -2*100*y(242) + 2*100*(24/50)*y(241); % plasma cysteine

dydt(243)= (Qliver/Vplasmatotal)*y(245) - (Qliver/Vplasmatotal)*y(243) - gsh_deg*y(243); %plasma gsh

dydt(244)= (Qliver/(Vcapliver*(1-hematocrit)))*y(242) - (Qliver/(Vcapliver*(1-hematocrit)))*y(244)  - fracliver*Vcysin(y(244)); %lcp cysteine  

dydt(245)= (Vtissliver*fracliver*0.85/(Vcapliver*(1-hematocrit)))*vGSHout_h(y(205)) + (Vtissliver*fracliver*0.85/(Vcapliver*(1-hematocrit)))*vGSHout_l(y(205)) - (Qliver/(Vcapliver*(1-hematocrit)))*y(245) + (Qliver/(Vcapliver*(1-hematocrit)))*y(243); %lcp gsh

dydt=real(dydt);

end
