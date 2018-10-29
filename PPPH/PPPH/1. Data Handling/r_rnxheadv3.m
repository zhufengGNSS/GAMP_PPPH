function [inf] = r_rnxheadv3(f_obs)

[fid,errmsg] = fopen(f_obs);

if any(errmsg)
    errordlg('OBSERVATION file can not be opened !','Observation File Error');
    error   ('OBSERVATION file can not be opened !');
end

inf.time.leap = [];

if size(f_obs,2)>28
    inf.time.int = str2double(f_obs(29:30));%sec
else
    inf.time.int  = 30;%sec
end

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
            inf.ant.number = strtrim(tline( 1:20));
            inf.ant.type   = strtrim(tline(21:40));
        case 'APPROX POSITION XYZ'
            
            inf.rec.pos    = sscanf(tline(1:60),'%f',[1,3]);
        case 'ANTENNA: DELTA H/E/N'
            inf.ant.hen = sscanf(tline(1:60),'%f',[1,3]);
        case 'SYS / # / OBS TYPES'
            if strcmp(tline(1),'G') %GPS
                no = sscanf(tline(5:6),'%d');
                inf.nob.gps = no;
                if no<14
                    lst = sscanf(tline(8:60),'%s');
                else
                    l1 = sscanf(tline(8:60),'%s');
                    tline = fgetl(fid);
                    l2 = sscanf(tline(8:60),'%s');
                    lst = strcat(l1,l2);
                end
                
                inf.seq.gps    = zeros(1,4);
                prior = 'PWCSLXYMND';
                % P1
                for pi = 1:10
                    if any(strfind(lst,strcat('C1',prior(pi))))
                        inf.seq.gps(1) = (strfind(lst,strcat('C1',prior(pi))) + 2)/3;
                        break
                    end
                end
                % P2
                for pi = 1:10
                    if any(strfind(lst,strcat('C2',prior(pi))))
                        inf.seq.gps(2) = (strfind(lst,strcat('C2',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L1
                for pi = 1:10
                    if any(strfind(lst,strcat('L1',prior(pi))))
                        inf.seq.gps(3) = (strfind(lst,strcat('L1',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L2
                for pi = 1:10
                    if any(strfind(lst,strcat('L2',prior(pi))))
                        inf.seq.gps(4) = (strfind(lst,strcat('L2',prior(pi))) + 2)/3;
                        break
                    end
                end
            elseif strcmp(tline(1),'R') %GLONASS
                no = sscanf(tline(5:6),'%d');
                inf.nob.glo = no;
                if no<14
                    lst = sscanf(tline(8:60),'%s');
                else
                    l1 = sscanf(tline(8:60),'%s');
                    tline = fgetl(fid);
                    l2 = sscanf(tline(8:60),'%s');
                    lst = strcat(l1,l2);
                end
                
                inf.seq.glo    = zeros(1,4);
                prior = 'PC';
                % P1
                for pi = 1:2
                    if any(strfind(lst,strcat('C1',prior(pi))))
                        inf.seq.glo(1) = (strfind(lst,strcat('C1',prior(pi))) + 2)/3;
                        break
                    end
                end
                % P2
                for pi = 1:2
                    if any(strfind(lst,strcat('C2',prior(pi))))
                        inf.seq.glo(2) = (strfind(lst,strcat('C2',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L1
                for pi = 1:2
                    if any(strfind(lst,strcat('L1',prior(pi))))
                        inf.seq.glo(3) = (strfind(lst,strcat('L1',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L2
                for pi = 1:2
                    if any(strfind(lst,strcat('L2',prior(pi))))
                        inf.seq.glo(4) = (strfind(lst,strcat('L2',prior(pi))) + 2)/3;
                        break
                    end
                end
            elseif strcmp(tline(1),'E') %GALILEO
                no = sscanf(tline(5:6),'%d');
                inf.nob.gal = no;
                if no<14
                    lst = sscanf(tline(8:60),'%s');
                else
                    l1 = sscanf(tline(8:60),'%s');
                    tline = fgetl(fid);
                    l2 = sscanf(tline(8:60),'%s');
                    lst = strcat(l1,l2);
                end
                
                inf.seq.gal    = zeros(1,4);
                prior1 = 'BCXAZ';
                prior5 = 'IQX';
                % P1
                for pi = 1:5
                    if any(strfind(lst,strcat('C1',prior1(pi))))
                        inf.seq.gal(1) = (strfind(lst,strcat('C1',prior1(pi))) + 2)/3;
                        break
                    end
                end
                % P5
                for pi = 1:3
                    if any(strfind(lst,strcat('C5',prior5(pi))))
                        inf.seq.gal(2) = (strfind(lst,strcat('C5',prior5(pi))) + 2)/3;
                        break
                    end
                end
                % L1
                for pi = 1:5
                    if any(strfind(lst,strcat('L1',prior1(pi))))
                        inf.seq.gal(3) = (strfind(lst,strcat('L1',prior1(pi))) + 2)/3;
                        break
                    end
                end
                % L5
                for pi = 1:3
                    if any(strfind(lst,strcat('L5',prior5(pi))))
                        inf.seq.gal(4) = (strfind(lst,strcat('L5',prior5(pi))) + 2)/3;
                        break
                    end
                end
            elseif strcmp(tline(1),'C') %BEIDOU
                no = sscanf(tline(5:6),'%d');
                inf.nob.bds = no;
                if no<14
                    lst = sscanf(tline(8:60),'%s');
                else
                    l1 = sscanf(tline(8:60),'%s');
                    tline = fgetl(fid);
                    l2 = sscanf(tline(8:60),'%s');
                    lst = strcat(l1,l2);
                end
                
                inf.seq.bds    = zeros(1,4);
                prior = 'IQX';
                % P1
                for pi = 1:3
                    if any(strfind(lst,strcat('C1',prior(pi))))
                        inf.seq.bds(1) = (strfind(lst,strcat('C1',prior(pi))) + 2)/3;
                        break
                    elseif any(strfind(lst,strcat('C2',prior(pi))))
                        inf.seq.bds(1) = (strfind(lst,strcat('C2',prior(pi))) + 2)/3;
                        break
                    end
                end
                % P7
                for pi = 1:3
                    if any(strfind(lst,strcat('C7',prior(pi))))
                        inf.seq.bds(2) = (strfind(lst,strcat('C7',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L1
                for pi = 1:3
                    if any(strfind(lst,strcat('L1',prior(pi))))
                        inf.seq.bds(3) = (strfind(lst,strcat('L1',prior(pi))) + 2)/3;
                        break
                    elseif any(strfind(lst,strcat('L2',prior(pi))))
                        inf.seq.bds(3) = (strfind(lst,strcat('L2',prior(pi))) + 2)/3;
                        break
                    end
                end
                % L7
                for pi = 1:3
                    if any(strfind(lst,strcat('L7',prior(pi))))
                        inf.seq.bds(4) = (strfind(lst,strcat('L7',prior(pi))) + 2)/3;
                        break
                    end
                end
            end
        case 'INTERVAL'
            inf.time.int = sscanf(tline(1:10),'%f');
        case 'TIME OF FIRST OBS'
            inf.time.first  = sscanf(tline( 1:44),'%d');
            inf.time.system = sscanf(tline(45:60),'%s');
        case 'TIME OF LAST OBS'
            inf.time.last   = sscanf(tline( 1:44),'%d');
        case 'LEAP SECONDS'
            inf.time.leap = sscanf(tline(1:6),'%d');
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