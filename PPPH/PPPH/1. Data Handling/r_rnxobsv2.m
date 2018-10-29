function [obs] = r_rnxobsv2(f_obs,inf,options,dcb)

[fid,errmsg] = fopen(f_obs);

if any(errmsg)
    errordlg('RINEX file can not be opened !','RINEX file error');
    error   ('RINEX file error');
end

sno = 105;

if inf.time.last == 86400
    fi = inf.time.first(4,1)*3600 + inf.time.first(5,1)*60 + inf.time.first(6,1);
    la = inf.time.last;
    max = (la-fi)/inf.time.int + 1;
else
    fi = inf.time.first(4,1)*3600 + inf.time.first(5,1)*60 + inf.time.first(6,1);
    la = inf.time.last (4,1)*3600 + inf.time.last (5,1)*60 + inf.time.last (6,1);
    max = (la-fi)/inf.time.int + 1;
end

p1s  =   NaN(max,sno);
p2s  =   NaN(max,sno);
l1s  =   NaN(max,sno);
l2s  =   NaN(max,sno);
eps  =   NaN(max, 1);
sts  =  zeros(max,sno);
[~,wavl] = frequencies;


d.year = inf.time.first(1);
d.mon  = inf.time.first(2);
d.day  = inf.time.first(3);
if d.year<2000
    d.year = d.year - 1900;
else
    d.year = d.year - 2000;
end


epno  = 0;
linenum = 0;
while ~feof(fid)
    
    tline = fgetl(fid);
    linenum = linenum + 1;
    
    if size(tline,2)<32
        ep    = sscanf(tline(1:end),'%f',[1,8]);
    else
        ep    = sscanf(tline(1:32),'%f',[1,8]);
    end
    
    if length(ep)==8 && ep(1)==d.year && ep(2)==d.mon && ep(3)==d.day && size(tline,2)>33
        epno  = epno + 1;
        epocs = ep(4)*3600 + ep(5)*60 + ep(6);
        eps(epno,1) = epocs;
        sats_no = ep(8);
        
        if sats_no<13
            if length(tline)<68
                satline = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline = strtrim(tline(33:68));
            end
        elseif sats_no<25
            if length(tline)<68
                satline1 = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline1 = strtrim(tline(33:68));
            end
            tline = fgetl(fid);
            linenum = linenum + 1;
            if length(tline)<68
                satline2 = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline2 = strtrim(tline(33:68));
            end
            satline = strcat(satline1,satline2);
        elseif sats_no<37
            if length(tline)<68
                satline1 = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline1 = strtrim(tline(33:68));
            end
            tline = fgetl(fid);
            linenum = linenum + 1;
            if length(tline)<68
                satline2 = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline2 = strtrim(tline(33:68));
            end
            tline = fgetl(fid);
            linenum = linenum + 1;
            if length(tline)<68
                satline3 = strtrim(tline(33:end));
            elseif length(tline)<=80
                satline3 = strtrim(tline(33:68));
            end
            satline = strcat(satline1,satline2,satline3);
        end
        
        if any(strfind(satline,'G')) || any(strfind(satline,'R')) || any(strfind(satline,'E'))
            satline=strrep(satline,'G','1');
            satline=strrep(satline,' ','0');
            satline=strrep(satline,'R','2');
            satline=strrep(satline,'E','3');
            satline=strrep(satline,'S','4');
            satline=strrep(satline,'T','4');
        end
        
        sat = sscanf(satline,'%3d');
        
        if options.system.gps == 1
            dtc = find(sat>100 & sat<133);
            sat(dtc) = sat(dtc) - 100;
        end
        
        if options.system.glo == 1
            dtc = find(sat>200 & sat<227);
            sat(dtc) = sat(dtc) - 200 + 32;
        end
        
        if options.system.gal == 1
            dtc = find(sat>310 & sat<331);
            sat(dtc) = sat(dtc) - 300 + 48;
        end
        
        for i=sat'
            
            ls = NaN(1,inf.obsno);
            if i>100
                for k=1:ceil(inf.obsno/5)
                    fgetl(fid);
                    linenum = linenum + 1;
                end
            else
                for k=1:ceil(inf.obsno/5)
                    tline = fgetl(fid);
                    linenum = linenum + 1;
                    for n=0:4
                        st = (16*n+1); fn = (16*(n+1)-2);
                        if fn<=length(tline)
                            tls= sscanf(tline(st:fn),'%f');
                            if ~isempty(tls)
                                ls((5*(k-1))+(n+1)) = tls;
                            end
                        end
                    end
                end
                if i<33
                    % P1
                    if (inf.seq.gps(1)==0) || (isnan(ls(inf.seq.gps(1))))
                        if options.dcb == 1
                            p1s(epno,i) = ls(inf.seq.gps(5)) + dcb(i,1);
                        else
                            p1s(epno,i) = ls(inf.seq.gps(5));
                        end
                    else
                        p1s(epno,i) = ls(inf.seq.gps(1));
                    end
                    % P2
                    p2s(epno,i) = ls(inf.seq.gps(2));
                    % L1
                    l1s(epno,i) = ls(inf.seq.gps(3))*wavl(i,1);
                    % L2
                    l2s(epno,i) = ls(inf.seq.gps(4))*wavl(i,2);
                elseif i<59
                    % P1
                    if (inf.seq.glo(1)==0) || (isnan(ls(inf.seq.glo(1))))
                        if options.dcb == 1
                            p1s(epno,i) = ls(inf.seq.glo(5)) + dcb(i,1);
                        else
                            p1s(epno,i) = ls(inf.seq.glo(5));
                        end
                    else
                        p1s(epno,i) = ls(inf.seq.glo(1));
                    end
                    p2s(epno,i) = ls(inf.seq.glo(2));
                    l1s(epno,i) = ls(inf.seq.glo(3))*wavl(i,1);
                    l2s(epno,i) = ls(inf.seq.glo(4))*wavl(i,2);
                elseif i<89
                    p1s(epno,i) = ls(inf.seq.gal(1));
                    p2s(epno,i) = ls(inf.seq.gal(2));
                    l1s(epno,i) = ls(inf.seq.gal(3))*wavl(i,1);
                    l2s(epno,i) = ls(inf.seq.gal(4))*wavl(i,2);
                end
            end
        end
    end
end

all = p1s + p2s + l1s + l2s;
sts(~isnan(all)) = 1;

if max>epno
    p1s(epno+1:max,:) = [];
    p2s(epno+1:max,:) = [];
    l1s(epno+1:max,:) = [];
    l2s(epno+1:max,:) = [];
    eps(epno+1:max,:) = [];
    sts(epno+1:max,:) = [];
end

obs.p1 = p1s;
obs.p2 = p2s;
obs.l1 = l1s;
obs.l2 = l2s;
obs.ep = eps;
obs.st = sts;
fclose('all');
end