% Rat glutathione synthetase (GS) (1)
% from Luo00 paper
% No feedback inhibition by GSH (Lu09)
% BRENDA does not list any kinetic parameters for mouse

function V = VGSSrat(glutamylcys);

atp=3000;
gly=1700;

km_gly=913;
km_atp=37;
km_glutamylcys=42;  %really for gluABA, which was used instead in Dinescu04
alpha=7.2;



vmax=24287;

V=vmax*(gly./(gly+km_gly)).*(atp./(atp+km_atp))...
.*(alpha*km_glutamylcys*glutamylcys+glutamylcys.^2)./(alpha*km_glutamylcys^2+2*alpha*km_glutamylcys*glutamylcys+glutamylcys.^2);








%k1 = 452; % 0.452 mM glycine Njalsson 2001
%k2 = 63; % 0.063 mM first binding glutamyl-cysteine (value from gamma-glutamyl-alpha-aminobutyrate) Njalsson 2001
%k3 = 164; % 0.164 mM second binding glutamyl-cysteine (value from gamma-glutamyl-alpha-aminobutyrate) Njalsson 2001
%V = (336)*3;  % 5.6umol/(min*mg enzyme) I think this should be high enough so that it will never be the rate limiting reaction

%V = 4500; %5000
%a = V.*((b.*c)./(k1.*k2) + (b.*c.*c)./(2.6.*k1.*k2.*k2))./(1 + (b./k1) + ((2.*c)./k2) + ((c.*c)/(2.6.*k2.*k2)) + ((2.*b.*c)./(k1.*k2)) +  ((b.*c.*c)/(2.6.*k1.*k2.*k2))); 
% test this eqn on monday
%--------------------------------------------------------------------------

% this enzyme shows negative cooperativity in the binding of
% glutamyl-cysteine, in that binding of the first substrate (k2) reduces enzyme
% affinity for the second substrate (k3) by a factor of alpha, which was
% calculated to be 2.6. The adair equation can be broken down into what
% resembles a two substrate M-M equation (1 substrate concentration but 2
% Kms), so I added in the "third" substrate, glycine, to make the equation
% work for three variables
