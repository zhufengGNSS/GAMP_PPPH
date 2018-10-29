function [atx] = r_antx(f_anten,inf,options)

[fid,errmsg] = fopen(f_anten);

if any(errmsg)
    errordlg('Antex file can not be opened.','Antex file error');
    error   ('Antex file error');
end

type = inf.ant.type;
sn   = 105;
sat_neu = NaN(sn,3,2);
rcv_neu = NaN( 1,3,4);

linenum = 0;
while ~feof(fid)
    
    tline = fgetl(fid);
    linenum = linenum + 1;
    tag   = strtrim(tline(61:end));
    
    if strcmp(tag,'START OF ANTENNA')
        tline = fgetl(fid);
        linenum = linenum + 1;
        tag   = strtrim(tline(61:end));
        
        if strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'G') && (options.system.gps == 1)
            sat_no = sscanf(tline(22:23),'%d');
            if sat_no>32
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'R') && (options.system.glo == 1)
            sat_no = 32 + sscanf(tline(22:23),'%d');
            if sat_no>58 
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'E') && (options.system.gal == 1)
            sat_no = 58 + sscanf(tline(22:23),'%d');
            if sat_no>88
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E05')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'C') && (options.system.bds == 1)
            sat_no = 88 + sscanf(tline(22:23),'%d');
            if sat_no>105
                continue
            end
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C07')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        sat_neu(sat_no,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end
            
        elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(strtrim(tline(1:20)),type)
            while ~strcmp(tag,'END OF ANTENNA')
                tline = fgetl(fid);
                linenum = linenum + 1;
                tag   = strtrim(tline(61:end));
                if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                    frq_no = 1; %L1
                    tline = fgetl(fid);
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                    frq_no = 2; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                    frq_no = 3; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                    frq_no = 4; %L2
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
                    if strcmp(tag,'NORTH / EAST / UP')
                        rcv_neu(1,:,frq_no) = sscanf(tline,'%f',[1,3]);
                    end
                end
            end 
        end
    else
        continue
    end
end

atx.sat.neu = sat_neu./1000;
atx.rcv.neu = rcv_neu./1000;

fclose('all');
end