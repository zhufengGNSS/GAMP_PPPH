function [Header,Epoch] = read_ppp(filename)
% SYNTAX:
%     =======================================
%     | [Header,Epoch] = read_ppp(filename) |
%     =======================================
%
%     read PPP intial file
%
% INPUTS:
% filename: PPP intial file
%
% OUTPUT:
%   Header: header information
%    Epoch: epoch data
%
% Originally written by Feng Zhou on 16/04/2016 @ GFZ
%
% Email: zhouforme@163.com
%
% Section 1.1, Space Geodetic Techniques, German Research Centre for Geosciences (GFZ)
%
%%===============================BEGIN PROGRAM=====================================%%
global DEG2RAD

fid = fopen(filename,'r');
if (fid == -1), error('*** ERROR(read_ppp): failed to open file: %s',filename); end

% read the header
on = 0;
while (1)
    line = fgets(fid);
    if (line == -1), break; end
    if (strcmp(deblank(line),'+PPP_HEADER'))
        on = 1;
        continue;
    end
    if (strcmp(deblank(line),'-PPP_HEADER')), break; end
    if (~on), continue; end
    switch deblank(line(1:16))
        case '       STA_NAME:'
            Header.stanam = deblank(line(18:end));   % station name
        case '       RCV_TYPE:'
            Header.rcvtyp = deblank(line(18:end));   % receiver type
        case '       ANT_TYPE:'
            Header.anttyp = deblank(line(18:end));   % antenna type
        case '        STA_POS:'
            Header.stapos = str2num(line(18:59));    % station position
        case '       INTERVAL:'
            Header.dintv = str2double(line(18:23));  % observation interval
        case '       GPS_TYPE:'
            Header.gpstyp = strtrim(line(18:end));   % GPS observation type
        case '   GLONASS_TYPE:'
            Header.glotyp = strtrim(line(18:end));   % GLONASS observation type
        case '    BEIDOU_TYPE:'
            Header.bdstyp = strtrim(line(18:end));   % BeiDou observation type
        case '   GALILEO_TYPE:'
            Header.galtyp = strtrim(line(18:end));   % Galileo observation type
        case '     TROP MODEL:'
            Header.trpmod = deblank(line(18:end));   % tropospheric model
        case '       TROP MAP:'
            Header.trpmap = deblank(line(18:end));   % tropospheric map function
        case '      ZTD MODEL:'
            Header.ztdmod = deblank(line(18:end));   % ZTD model
        case '     ORBIT TYPE:'
            Header.orbtyp = deblank(line(18:end));   % orbit type
        case '     CLOCK TYPE:'
            Header.clktyp = deblank(line(18:end));   % clock type
        case '     ELEV. MASK:'
            Header.elemask = str2double(line(18:21))*DEG2RAD;  % elevation mask
        case '      GLONASS K:'
            Header.glofrq = str2num(line(18:end));   % GLONASS frequency number
    end
end

% read the block
on = 0;
iepo = 0;
while (1)
    line = fgets(fid);
    if (line == -1), break; end
    if (strcmp(deblank(line),'+PPP_BLOCK'))
        on = 1;
        continue;
    end
    if (strcmp(deblank(line),'-PPP_BLOCK')), break; end
    if (~on), continue; end
    if (line(1:1) ~= '>'), continue; end
    iepo = iepo+1;
    Temp = str2num(line(2:end));
    if (length(Temp) ~= 7), error(' *** ERROR(read_ppp): error in epoch line'); end
    Epoch(iepo).nsat = Temp(7);   % number of sat.
    Epoch(iepo).tepo.jd = transT_ymd2mjd(Temp(1),Temp(2),Temp(3));   % epoch time tag
    Epoch(iepo).tepo.sod = transT_hms2sod(Temp(4),Temp(5),Temp(6));
    for ii = 1:Epoch(iepo).nsat
        line = fgets(fid);
        Epoch(iepo).sats(ii,1) = transS_satid2no(line(2:4));   % sat. number
        Epoch(iepo).csats{ii,1} = line(2:4);   % sat. id
        Epoch(iepo).flag(ii,1) = str2double(line(6:6));   % flag
        Epoch(iepo).lam(ii,1:3) = get_satwavelen(Epoch(iepo).sats(ii),Header.glofrq);   % sat. carrier wave length
        Epoch(iepo).satPos(ii,1:3) = str2num(line(8:54));     % sat. position (m)
        Epoch(iepo).satClk(ii,1) = str2double(line(56:70));   % sat. clock (m)
        Epoch(iepo).elev(ii,1) = str2double(line(72:86))*DEG2RAD;     % sat. elevation (rad)
        Epoch(iepo).azim(ii,1) = str2double(line(88:102))*DEG2RAD;    % sat. azimuth (rad)
        Epoch(iepo).P(ii,1) = str2double(line(104:118));      % P1 code observation (m)
        Epoch(iepo).P(ii,2) = str2double(line(120:134));      % P2 code observation (m)
        Epoch(iepo).P(ii,3) = str2double(line(136:150));      % P2 code observation (m)
        Epoch(iepo).L(ii,1) = str2double(line(152:166));      % L1 phase observation (cycles)
        Epoch(iepo).L(ii,2) = str2double(line(168:182));      % L2 phase observation (cycles)
        Epoch(iepo).L(ii,3) = str2double(line(184:198));      % L3 phase observation (cycles)
        Epoch(iepo).dtrp(ii,1) = str2double(line(200:214));   % slant tropospheric delay (m)
        Epoch(iepo).wmap(ii,1) = str2double(line(216:230));   % slant tropospheric wet map function (m)
        Epoch(iepo).dsag(ii,1) = str2double(line(232:246));   % sagnac effet (m)
        Epoch(iepo).dtid(ii,1) = str2double(line(248:262));   % tide effect (solid earth tide+pole tide+ocean tide, m)
        Epoch(iepo).pcvl1(ii,1) = str2double(line(264:278));  % PCV on L1 (m)
        Epoch(iepo).pcvl2(ii,1) = str2double(line(280:294));  % PCV on L2 (m)
        Epoch(iepo).windup(ii,1) = str2double(line(296:310)); % phase windup (cycles)
    end
end

return;
%
%%================================END PROGRAM======================================%%