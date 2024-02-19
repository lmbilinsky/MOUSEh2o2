% we assume VGR doesn't depend on GSH

function a = VGRmito(gssg_mito, cNADPH, gs_apap);

km1 = 107;        % for gssg
km2 = 10.4;       % for NADPH
km3 = 3000;       % gs_apap




vmax=42948;






a = vmax*((gssg_mito/km1).*(cNADPH/km2))./(1 + gssg_mito/km1 + cNADPH/km2 + (gssg_mito/km1).*(cNADPH/km2));

