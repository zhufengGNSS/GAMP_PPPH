function lam = get_satwavelen(sat,glo_frq)
% SYNTAX:
%     =====================================
%     | lam = get_satwavelen(sat,glo_frq) |
%     =====================================
%
%     To get GNSS satellite carrier wave lengths
%
% INPUTS:
%       sat: satellite number
%   glo_frq: GLONASS frequency numbers
%
% OUTPUT:
%       lam: carrier wave lengths (m)
%
% Originally written by Feng Zhou on 16/04/2016 @ GFZ
%
% Email: fzhou@geodesy.cn; fzhou@gfz-potsdam.de; zhouforme@gmail.com
%
% Section 1.1, Space Geodetic Techniques, German Research Centre for Geosciences (GFZ)
%
%%===============================BEGIN PROGRAM=====================================%%
global FREQ1 FREQ2 FREQ5 FREQ7
global FREQ1_GLO DFRQ1_GLO FREQ2_GLO DFRQ2_GLO FREQ3_GLO
global FREQ1_BDS FREQ2_BDS FREQ3_BDS
global SYS_GPS SYS_GLO SYS_BDS SYS_GAL
global CLIGHT

[sys,prn] = transS_satsys(sat);
if (sys == SYS_GPS)      % GPS
    lam1 = CLIGHT/FREQ1;
    lam2 = CLIGHT/FREQ2;
    lam3 = CLIGHT/FREQ5;
elseif (sys == SYS_GLO)  % GLONASS
    lam1 = CLIGHT/(FREQ1_GLO+DFRQ1_GLO*glo_frq(prn));
    lam2 = CLIGHT/(FREQ2_GLO+DFRQ2_GLO*glo_frq(prn));
    lam3 = CLIGHT/FREQ3_GLO;
elseif (sys == SYS_BDS)  % BeiDou
    lam1 = CLIGHT/FREQ1_BDS;
    lam2 = CLIGHT/FREQ2_BDS;
    lam3 = CLIGHT/FREQ3_BDS;
elseif (sys == SYS_GAL)  % Galileo
    lam1 = CLIGHT/FREQ1;
    lam2 = CLIGHT/FREQ5;
    lam3 = CLIGHT/FREQ7;
end
lam = [lam1 lam2 lam3];

return;
%
%%================================END PROGRAM======================================%%