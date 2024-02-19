function a = VGPXcyt(gsh_cyt, H2O2);

% we assume VGPX doesn't depend on GSSG
km_g = 1330; %see page 53 of Carmagnol83. Very close to the value on page 410 of Yamashita12    
km_h = 6.8; %from Yamashita12 (GPx1 in red muscle of bluefin tuna. GPx1 has preferred subsrate H2O2)

vmax=4332453;

a = vmax .* ((gsh_cyt./(km_g + gsh_cyt))^2) .* (H2O2./(km_h + H2O2));

