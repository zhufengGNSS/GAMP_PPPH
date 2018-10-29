function [sat,inf] = r_sp3(f_orbit,options,inf)

[fid,errmsg] = fopen(f_orbit);

if any(errmsg)
    errordlg('SP3 file can not be opened.','SP3 file error');
    error   ('SP3 file error');
end

sn  = 105;
sp3 = NaN(96,4,sn); 

while ~feof(fid)
    line = fgetl(fid);
    if strcmp(line(1),'#') && ~strcmp(line(1:2),'##')
        dat = sscanf(line(4:14),'%f');
    end
    if strcmp(line(1:2),'##')
        sp3int = sscanf(line(25:38),'%f');
        inf.time.sp3int = sp3int;
        epno   = 86400/sp3int; 
        sp3    = NaN(epno,4,sn);
    end
    
    if line(1)=='+'
        temp = sscanf(line(5:6),'%d');
        if ~isnan(temp)
            NoSat = temp;
        end
    end
    
    % new epoch
    if line(1)=='*'
        ep   = sscanf(line(2:end),'%f',[1,6]);
        if ep(1)~=dat(1) || ep(2)~=dat(2) || ep(3)~=dat(3)
            continue
        end
        epno = ((ep(4)*3600 + ep(5)*60 + ep(6))/sp3int)+1;
        for k=1:NoSat
            line = fgetl(fid);
            if strcmp(line(2),'G') && (options.system.gps == 1)
                sno  = sscanf(line(3:4),'%d');
                if sno>32
                    continue
                else
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'R') && (options.system.glo == 1)
                sno  = 32 + sscanf(line(3:4),'%d');
                if sno>58
                    continue
                else
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'E') && (options.system.gal == 1)
                sno  = 58 + sscanf(line(3:4),'%d');
                if sno>88
                    continue
                else
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    % writing part
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            elseif strcmp(line(2),'C') && (options.system.bds == 1)
                sno  = 88 + sscanf(line(3:4),'%d');
                if sno>105
                    continue
                else
                    temp = sscanf(line(5:end),'%f',[1,4]);
                    
                    sp3(epno,1:3,sno) = temp(1:3)*1000; %meter
                    sp3(epno,  4,sno) = temp(4)*10^-6;  %second 
                end
            end        
        end        
    end
end

sat.sp3    = sp3;

fclose('all');
end
