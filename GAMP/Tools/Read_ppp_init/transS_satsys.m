function [sys,prn] = transS_satsys(sat)
% SYNTAX:
%     ==================================
%     | [sys,prn] = transS_satsys(sat) |
%     ==================================
%
%     convert satellite number to satellite system
%
% INPUTS:
%     sat: satellite number
%
% OUTPUT:
%     sys: satellite system
%     prn: satellite prn
%
% Originally written by Feng Zhou on 18/04/2016 @ GFZ
%
% Email: fzhou@geodesy.cn; fzhou@gfz-potsdam.de; zhouforme@gmail.com
%
% Section 1.1, Space Geodetic Techniques, German Research Centre for Geosciences (GFZ)
%
%%===============================BEGIN PROGRAM=====================================%%
global SYS_NONE SYS_GPS SYS_GLO SYS_BDS SYS_GAL
global NSAT_GPS NSAT_GLO NSAT_BDS NSAT_GAL MAXSAT

sys = SYS_NONE;
if (sat <= 0 || sat > MAXSAT)
    sat = 0;
elseif (sat <= NSAT_GPS)
    sys = SYS_GPS;
elseif (sat-NSAT_GPS <= NSAT_GLO)
    sys = SYS_GLO;
    sat = sat-NSAT_GPS;
elseif (sat-NSAT_GPS-NSAT_GLO <= NSAT_BDS)
    sys = SYS_BDS;
    sat = sat-NSAT_GPS-NSAT_GLO;
elseif (sat-NSAT_GPS-NSAT_GLO-NSAT_BDS <= NSAT_GAL)
    sys = SYS_GAL;
    sat = sat-NSAT_GPS-NSAT_GLO-NSAT_BDS;
else
    sat = 0;
end
prn = sat;

return;
%
%%================================END PROGRAM======================================%%