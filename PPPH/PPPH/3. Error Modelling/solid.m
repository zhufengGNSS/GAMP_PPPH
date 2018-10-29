function [stide] = solid(r_xyz,s_xyz,sun_xyz,mon_xyz)

l   = r_xyz - s_xyz;
los = l./norm(l);

h0 = 0.6078; h2 = -0.0006; h3 = 0.292;
l0 = 0.0847; l2 =  0.0002; l3 = 0.015;

MS2E = 332946.0; MM2E = 0.01230002; re = 6378137;

[elip] = xyz2plh(r_xyz,1);
lat = elip(1);
trm = 3*sind(lat)^2 - 1;
h = h0 + h2*trm;
l = l0 + l2*trm;

sunDist  = norm(sun_xyz);
moonDist = norm(mon_xyz);
recDist  = norm(r_xyz);

sunUni  = sun_xyz./sunDist;
moonUni = mon_xyz./moonDist; 
recUni  = r_xyz./recDist;

dotSR = dot(sunUni,recUni);
dotMR = dot(moonUni,recUni);

aSun  =  sunUni - (dotSR.*recUni);
aMoon = moonUni - (dotMR.*recUni);

DRS  = (re^4)/(sunDist^3);
DRM  = (re^4)/(moonDist^3);
DRM2 = (re^5)/(moonDist^4);

s1  = ((3*(dotSR^2)) - 1)/2;
d2s = (MS2E*DRS).*(((h.*recUni).*s1) + ((3*l*dotSR).*aSun ));
s2  = ((3*(dotMR^2)) - 1)/2;
d2m = (MM2E*DRM).*(((h.*recUni).*s2) + ((3*l*dotMR).*aMoon));

s3  = ((5/2)*(dotMR^3)) - ((3/2)*(dotMR));
s4  = ((15/2)*(dotMR^2))- (3/2);
d3m = (MM2E*DRM2).*(((h3.*recUni).*s3) + ((l3*s4).*aMoon));

stide = d2s + d2m + d3m ;
stide = dot(stide,los);

end

