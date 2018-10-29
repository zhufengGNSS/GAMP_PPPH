function sat = transS_satid2no(csat)
% SYNTAX:
%     ===============================
%     | sat = transS_satid2no(csat) |
%     ===============================
%
%     convert satellite id to satellite number
%
% INPUTS:
%    csat: satellite id (Gnn,Rnn,Cnn,Enn)
%
% OUTPUT:
%     sat: satellite number
%
% Originally written by Feng Zhou on 18/04/2016 @ GFZ
%
% Email: fzhou@geodesy.cn; fzhou@gfz-potsdam.de; zhouforme@gmail.com
%
% Section 1.1, Space Geodetic Techniques, German Research Centre for Geosciences (GFZ)
%
%%===============================BEGIN PROGRAM=====================================%%
global NSAT_GPS NSAT_GLO NSAT_BDS

if (csat(1:1) == 'G')
    sat = str2double(csat(2:3));
elseif (csat(1:1) == 'R')
    sat = str2double(csat(2:3))+NSAT_GPS;
elseif (csat(1:1) == 'C')
    sat = str2double(csat(2:3))+NSAT_GPS+NSAT_GLO;
elseif (csat(1:1) == 'E')
    sat = str2double(csat(2:3))+NSAT_GPS+NSAT_GLO+NSAT_BDS;
end

return;
%
%%================================END PROGRAM======================================%%