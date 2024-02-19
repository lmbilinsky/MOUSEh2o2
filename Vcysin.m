%transport of cysteine into the cell from blood
%b = blood cysteine

function a = Vcysin(b);




Vmax=265132600;

%Vmax=0.0009*Vmax; %for arsenic model

km = 2100; %ASC transporter (Kilberg80)


a = Vmax*b./(km + b);



