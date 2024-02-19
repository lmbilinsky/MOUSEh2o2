% we assume VGR doesn't depend on GSH

function a = VGRcyt(gssg_cyt, cNADPH, gs_apap);

km1 = 107;        % for gssg
km2 = 10.4;       % for NADPH
km3 = 3000;       % gs_apap




vmax=12496;






a = vmax*((gssg_cyt/km1).*(cNADPH/km2))./(1 + gssg_cyt/km1 + cNADPH/km2 + (gssg_cyt/km1).*(cNADPH/km2));

