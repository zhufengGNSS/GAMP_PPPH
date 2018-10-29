function [dcb] = r_dcb(f_dcb)

[fid,errmsg] = fopen(f_dcb);

if any(errmsg)
    errordlg('DCB file can not be opened !','DCB File Error');
end


dcb = zeros(58,1);

while ~feof(fid)
    tline = fgetl(fid);
    if ~isempty(tline) && strcmp(tline(1),'G') && ~isempty(sscanf(tline(2:3),'%d'))
        li = sscanf(tline(2:end),'%f');
        if ~isempty(li)
            sno = li(1);
            if sno<33
                dcb(sno) = li(2);
            end
        else
            continue
        end
    end
    
    if ~isempty(tline) && strcmp(tline(1),'R') && ~isempty(sscanf(tline(2:3),'%d'))
        li = sscanf(tline(2:end),'%f');
        if ~isempty(li)
            sno = 32 + li(1);
            if sno<59
                dcb(sno) = li(2);
            end
        else
            continue
        end
    end
end

c = 299792458; %m/s
dcb = dcb.*(10^-9*c);

fclose('all');
end

