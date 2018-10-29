% SYNTAX:
%     ==================================
%     |  |
%     ==================================
%
%     to set global constant variables
%
% INPUTS:
%     
%
% OUTPUT:
%       
%
% Originally written by Feng Zhou on 16/04/2016 @ GFZ
%
% Email: fzhou@geodesy.cn; fzhou@gfz-potsdam.de; zhouforme@gmail.com
%
% Section 1.1, Space Geodetic Techniques, German Research Centre for Geosciences (GFZ)
%
%%===============================BEGIN PROGRAM=====================================%%
global FREQ1 FREQ2 FREQ5 FREQ6 FREQ7 FREQ8
global FREQ1_GLO DFRQ1_GLO FREQ2_GLO DFRQ2_GLO FREQ3_GLO
global FREQ1_BDS FREQ2_BDS FREQ3_BDS
global NSAT_GPS NSAT_GLO NSAT_BDS NSAT_GAL MAXSAT
global SYS_NONE SYS_GPS SYS_GLO SYS_BDS SYS_GAL
global OMGE_GPS OMGE_GLO OMGE_BDS OMGE_GAL
global CLIGHT DEG2RAD RAD2DEG


% GNSS satellite frequencies
FREQ1 = 1.57542e+9;       % GPS L1 or Galileo E1 frequency (Hz)
FREQ2 = 1.22760e+9;       % GPS L2 frequency (Hz)
FREQ5 = 1.17645e+9;       % GPS L5 or Galileo E5a frequency (Hz)
FREQ6 = 1.27875e+9;       % Galileo E6 frequency (Hz)
FREQ7 = 1.20714e+9;       % Galileo E5b frequency (Hz)
FREQ8 = 1.191795e+9;      % Galileo E5a+b frequency (Hz)
FREQ1_GLO = 1.60200e+9;   % GLONASS G1 base frequency (Hz)
DFRQ1_GLO = 0.56250e+6;   % GLONASS G1 bias frequency (Hz/n)
FREQ2_GLO = 1.24600e+9;   % GLONASS G2 base frequency (Hz)
DFRQ2_GLO = 0.43750e+6;   % GLONASS G2 bias frequency (Hz/n)
FREQ3_GLO = 1.202025e+9;  % GLONASS G3 frequency (Hz)
FREQ1_BDS = 1.561098e+9;  % BDS B1 frequency (Hz)
FREQ2_BDS = 1.20714e+9;   % BDS B2 frequency (Hz)
FREQ3_BDS = 1.26852e+9;   % BDS B3 frequency (Hz)

% satellite number
NSAT_GPS = 32;
NSAT_GLO = 27;
NSAT_BDS = 35;
NSAT_GAL = 30;
MAXSAT = NSAT_GPS+NSAT_GLO+NSAT_BDS+NSAT_GAL;

% satellite system
SYS_NONE = 0;
SYS_GPS  = 1;
SYS_GLO  = 2;
SYS_BDS  = 4;
SYS_GAL  = 8;

% other constants
OMGE_GPS = 7.2921151467e-5;   % earth angular velocity (IS-GPS) (rad/s)
OMGE_GLO = 7.292115e-5;       % earth angular velocity (rad/s)
OMGE_BDS = 7.292115e-5;       % earth angular velocity (rad/s)
OMGE_GAL = 7.2921151467e-5;   % earth angular velocity (rad/s)
CLIGHT = 299792458.0;         % speed of light (m/s)
DEG2RAD = pi/180.0;           % deg to rad
RAD2DEG = 180.0/pi;           % rad to deg


return;
%
%%================================END PROGRAM======================================%%