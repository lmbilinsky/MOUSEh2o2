function PrintFinalValuesMousePBPK(Vplasmatotal, Vtissliver, Vcapliverplasma, M, cysinput)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT OUT CONCENTRATIONS AND VELOCITIES AT FINAL TIME
t=datestr(now); % today's date
fid=fopen('FinalState.txt','w');  %this opens a file and then the vector components are put in  it.
fprintf(fid, '%s \n\n', t);



fprintf(fid, 'CON.(cell) \t\t\t VEL.(cell) \t\t\t CON.(lcp) \t\t VEL.(lcp) \t\t\t\t CON.(plasma) \t\t\t VEL.(plasma) \n\n');

  

fprintf(fid, 'H2O2cyt= %4.3f \t\t VGPXcyt=%4.2f \n',M(99),M(30));
fprintf(fid,' \t\t\t\t\t VCAT=%4.2f \n\n',M(113));
H2O2out=M(30)+M(113);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f \n\n',-H2O2out);
fprintf(fid, ' \t\t\t\t\t vh2o2other=%4.2f     \n',M(103));
fprintf(fid, ' \t\t\t\t\t diffh2o2mito_to_cyt=%4.2f     \n',  (0.15/0.85)*M(98));
H2O2in=(0.15/0.85)*M(98)+M(103);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f \n\n',H2O2in);                                   
                           
                              
fprintf(fid, 'H2O2mito= %4.3f \t VGPXmito=%4.2f \n',M(91),M(53));
fprintf(fid, ' \t\t\t\t\t diffh2o2mito_to_cyt=%4.2f \n',M(98));
H2O2out=M(53)+M(98);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f \n\n',-H2O2out);
fprintf(fid, ' \t\t\t\t\t VH2O2prod=%4.2f     \n',  M(89));
H2O2in=M(89);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f \n\n',H2O2in);                                   




fprintf(fid, 'lCYS= %4.2f \t\t VCYSDIOXY= %4.2f \t\t lcpCYS= %4.2f \t vlcpCYSl= %4.2f \t\t plasmaCYS= %4.2f \t\t vplasmaCYSlcp= %4.2f \n',M(4),M(17),M(111),M(27),M(109),M(106));
fprintf(fid, ' \t\t\t\t\t VGCL= %4.2f \t\t\t\t\t\t\t vlcpCYSplasma= %4.2f                     \n',M(19),M(104));
CYSoutliver=M(17)+M(19);
CYSoutlcp=M(104)+M(27);
CYSoutplasma=M(106);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f  \t\t\t\t\t\t TOTout=%4.2f \t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-CYSoutliver,-CYSoutlcp,-CYSoutplasma);
fprintf(fid, ' \t\t\t\t\t vlcpCYSl=%4.2f  \t\t\t\t\t\t vplasmaCYSlcp=%4.2f \t\t\t\t\t\t\t VCYSe=%4.2f    \n',(Vcapliverplasma/Vtissliver)*M(27),(Vplasmatotal/Vcapliverplasma)*M(106),cysinput);
fprintf(fid, ' \t\t\t\t\t VCYSmet=%4.2f \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t         VbGSHd=%4.2f                     \n',100,M(33));
fprintf(fid, ' \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t vlcpCYSplasma=%4.2f     \n',(Vcapliverplasma/Vplasmatotal)*M(104));
CYSinliver=(Vcapliverplasma/Vtissliver)*M(27)+100;
CYSinlcp=(Vplasmatotal/Vcapliverplasma)*M(106);
CYSinplasma=cysinput+M(33)+(Vcapliverplasma/Vplasmatotal)*M(104);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t\t\t\t TOTin=%4.2f \n\n',CYSinliver,CYSinlcp,CYSinplasma);


fprintf(fid, 'GluCys=%4.2f \t\t VGSS=%4.2f \n',M(6),M(20));
GluCysout=M(20);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f  \n\n',-GluCysout);
fprintf(fid, '\t\t\t\t     VGCL=%4.2f  \n',M(19)/0.85);
GluCysin=M(19)/0.85;
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f   \n\n\n',GluCysin);


fprintf(fid, 'cytGSH= %4.2f \t gshcyttomito=%4.2f   lcpGSH=%4.2f  \t vlcpGSHplasma=%4.2f   plasmaGSH=%4.2f \t\t vplasmaGSHlcp=%4.2f  \n',M(5),M(51),M(112),M(105),M(110),M(107));
fprintf(fid, ' \t\t\t\t\t VlGSHlcp=%4.2f \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t VbGSHd=%4.2f \n',M(28),M(33));
fprintf(fid, ' \t\t\t\t\t cytVGPX=%4.2f  \n',2*M(30));
GSHoutliver=M(51)+M(28)+2*M(30);
GSHoutlcp=M(105);
GSHoutplasma=M(107)+M(33);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f \t\t\t\t\t\t TOTout=%4.2f  \t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-GSHoutliver,-GSHoutlcp,-GSHoutplasma);

fprintf(fid, ' \t\t\t\t\t VGSS=%4.2f \t\t\t\t\t\t\t VlGSHlcp=%4.2f \t\t\t\t\t\t\t\t vlcpGSHplasma=%4.2f \n\n', M(20),(Vtissliver/Vcapliverplasma)*0.85*M(28),(Vcapliverplasma/Vplasmatotal)*M(105));
fprintf(fid, ' \t\t\t\t\t cytVGR= %4.2f \t\t\t\t\t\t vplasmaGSHlcp=%4.2f  \n',2*M(31),(Vplasmatotal/Vcapliverplasma)*M(107)); 
fprintf(fid, ' \t\t\t\t\t gshmitotocyt= %4.2f   \n',M(55)*(0.15/0.85));
GSHinliver=M(20)+2*M(31)+M(55)*(0.15/0.85);
GSHinlcp=(Vtissliver/Vcapliverplasma)*0.85*M(28)+(Vplasmatotal/Vcapliverplasma)*M(107);
GSHinplasma=(Vcapliverplasma/Vplasmatotal)*M(105);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t\t\t TOTin=%4.2f \n\n',GSHinliver,GSHinlcp,GSHinplasma);


fprintf(fid, 'mitoGSH= %4.2f  \t mitoVGPX=%4.2f \n',M(57),2*M(53));
fprintf(fid, ' \t\t\t\t\t gshmitotocyt= %4.2f   \n',M(55));
mitoGSHoutliver=2*M(53)+M(55);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f \n\n',-mitoGSHoutliver);
fprintf(fid, ' \t\t\t\t\t mitoVGR= %4.2f    \n',2*M(54));
fprintf(fid, ' \t\t\t\t\t gshcyttomito= %4.2f   \n',M(51)*(0.85/0.15));
mitoGSHinliver=2*M(54)+M(51)*(0.85/0.15);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f  \n\n\n',mitoGSHinliver);

fprintf(fid, 'cytGSSG=%4.2f \t\t cytVGR=%4.2f \n',M(12),M(31));
fprintf(fid, '\t\t\t\t     VlGSSGb=%4.2f  \n',M(29));
GSSGoutliver=M(31)+M(29);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f  \n\n',-GSSGoutliver);
fprintf(fid, '\t\t\t\t\t cytVGPX=%4.2f\n', M(30));
GSSGinliver=M(30);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f   \n\n\n',GSSGinliver);

fprintf(fid, 'mitoGSSG=%4.2f \t\t mitoVGR=%4.2f \n',M(58),M(54));
mitoGSSGoutliver=M(54);
fprintf(fid, ' \t\t\t\t\t TOTout=%4.2f  \n\n',-mitoGSSGoutliver);
fprintf(fid, '\t\t\t\t\t mitoVGPX=%4.2f\n', M(53));
mitoGSSGinliver=M(53);
fprintf(fid, ' \t\t\t\t\t TOTin=%4.2f   \n\n\n',mitoGSSGinliver);





fclose(fid);
open FinalState.txt;   %this opens a window with the info in it
  
