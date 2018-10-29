function [inf] = r_rnxheadv2(f_obs)

[fid,errmsg] = fopen(f_obs);

if any(errmsg)
    errordlg('OBSERVATION file can not be opened !','Observation File Error');
    error   ('OBSERVATION file can not be opened !');
end

inf.time.leap = [];

inf.time.int  = 30; %sec

inf.time.last = 86400; %sec

while 1
    
    tline = fgetl(fid);
    tag  = strtrim(tline(61:end));
    switch tag
        case 'RINEX VERSION / TYPE'
            
            if strcmp(sscanf(tline(21),'%c'),'O')
                inf.rinex.type = sscanf(tline(21),'%c');
            else
                errordlg('It is not a observation file !','Observation file error');
                error   ('It is not a observation file !');
            end
            
            inf.sat.system = sscanf(tline(41),'%c');
        case 'REC # / TYPE / VERS'
            
            inf.rec.number  = strtrim(tline( 1:20));
            inf.rec.type    = strtrim(tline(21:40));
            inf.rec.version = strtrim(tline(41:60));
        case 'ANT # / TYPE'
            
            inf.ant.number  = strtrim(tline( 1:20));
            inf.ant.type    = strtrim(tline(21:40));
        case 'APPROX POSITION XYZ'
            
            inf.rec.pos     = sscanf(tline,'%f',[1,3]);
        case 'ANTENNA: DELTA H/E/N'
            
            inf.ant.hen     = sscanf(tline,'%f',[1,3]);
        case '# / TYPES OF OBSERV'
            obno = sscanf(tline(1:6),'%d');
            if obno<=9
                temp = sscanf(tline(7:60),'%s');
            elseif obno<=18
                temp = sscanf(tline(7:60),'%s');
                tline = fgetl(fid);
                temp  = strcat(temp,sscanf(tline(7:60),'%s'));
            elseif obno<=27
                temp = sscanf(tline(7:60),'%s');
                tline = fgetl(fid);
                temp  = strcat(temp,sscanf(tline(7:60),'%s'));
                tline = fgetl(fid);
                temp  = strcat(temp,sscanf(tline(7:60),'%s'));
            end
            inf.obsno = obno;
            
            inf.seq.gps = zeros(1,5);
            inf.seq.glo = zeros(1,5);
            inf.seq.gal = zeros(1,4);
            % C1
            if ~isempty(strfind(temp,'C1'))
                inf.seq.gps(5) = (strfind(temp,'C1')+1)/2;
                inf.seq.glo(5) = (strfind(temp,'C1')+1)/2;
                inf.seq.gal(1) = (strfind(temp,'C1')+1)/2;
            end
            % P1
            if ~isempty(strfind(temp,'P1'))
                inf.seq.gps(1) = (strfind(temp,'P1')+1)/2;
                inf.seq.glo(1) = (strfind(temp,'P1')+1)/2;
            end
            % P2
            if ~isempty(strfind(temp,'P2'))
                inf.seq.gps(2) = (strfind(temp,'P2')+1)/2;
                inf.seq.glo(2) = (strfind(temp,'P2')+1)/2;
            end
            % C5
            if ~isempty(strfind(temp,'C5'))
                inf.seq.gal(2) = (strfind(temp,'C5')+1)/2;
            end
            % L1
            if ~isempty(strfind(temp,'L1'))
                inf.seq.gps(3) = (strfind(temp,'L1')+1)/2;
                inf.seq.glo(3) = (strfind(temp,'L1')+1)/2;
                inf.seq.gal(3) = (strfind(temp,'L1')+1)/2;
            end
            % L2
            if ~isempty(strfind(temp,'L2'))
                inf.seq.gps(4) = (strfind(temp,'L2')+1)/2;
                inf.seq.glo(4) = (strfind(temp,'L2')+1)/2;
            end
            % L5
            if ~isempty(strfind(temp,'L5'))
                inf.seq.gal(4) = (strfind(temp,'L5')+1)/2;
            end
        case 'INTERVAL'
            inf.time.int = sscanf(tline(1:10),'%d');
        case 'TIME OF FIRST OBS'
            inf.time.first  = sscanf(tline(1:43),'%d');
            inf.time.system = strtrim(tline(44:60));
        case 'TIME OF LAST OBS'
            inf.time.last   = sscanf(tline(1:43),'%d');
        case 'LEAP SECONDS'
            inf.time.leap = sscanf(tline,'%d');
        case 'END OF HEADER'
            break
    end
end

if isempty(inf.time.leap)
    [~,mjd] = cal2jul(inf.time.first(1),inf.time.first(2),inf.time.first(3),...
        (inf.time.first(4)*3600 + inf.time.first(5)*60 + inf.time.first(6)));
    
    inf.time.leap = leapsec(mjd);
    if strcmp(inf.time.system,'GPS')
        inf.time.leap = inf.time.leap - 19;
    end
end


[doy] = clc_doy(inf.time.first(1),inf.time.first(2),inf.time.first(3));
inf.time.doy = doy;

fclose('all');
end