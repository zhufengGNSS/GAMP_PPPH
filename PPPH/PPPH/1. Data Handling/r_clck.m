function [clk] = r_clck(f_clock,options)

[fid,errmsg] = fopen(f_clock);

if any(errmsg)
    errordlg('Clock file can not be opened.','Clock file error');
    error   ('Clock file error');
end

sn = 105;
tn = 86400/options.clck_int;
clk = NaN(tn,sn);

while ~feof(fid)
    tline = fgetl(fid);
    
    if strcmp(tline(1:4),'AS G') && (options.system.gps == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = new(1);
        if sat_no<33
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS R') && (options.system.glo == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 32 + new(1);
        if sat_no<59
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS E') && (options.system.gal == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 58 + new(1);
        if sat_no<89
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    elseif strcmp(tline(1:4),'AS C') && (options.system.bds == 1)
        new = sscanf(tline(5:end),'%f',[1,10]);
        sat_no = 88 + new(1);
        if sat_no<106
            epoch  = (new(5)*3600 + new(6)*60 + new(7))/options.clck_int + 1;
            clk(epoch,sat_no) = new(9);
        else
            continue
        end
    end 
end
fclose('all');
end