function [cvg_n,cvg_e,cvg_u,rms_n,rms_e,rms_u] = analysis(pathname,filename,mode)

% NCOVER = 60;
NCOVER = 20;

% if filename(5) == '_'
%     staname = sprintf('%s%s%s%s%s%s',filename(1:4),' + ',filename(6:13), ' + ',mode,' + sp3');
% else
%     staname = sprintf('%s%s%s%s%s%s',filename(1:4),' + ',filename(7:14), ' + ',mode,' + realtime');
% end;
fullname = [pathname  filename];
A=load(fullname);
TIME = A(:,7:8);
ENU = A(:,12:14);
[m,n]=size(ENU);
success = 0;
nepo = 0;
cvg_n = m;
for i=1:m - NCOVER
    if abs(ENU(i,2)) < 0.1
        for j=i:i+NCOVER -1
            if abs(ENU(j,2)) < 0.1
                nepo = nepo + 1;
            end;
        end;
        if nepo == NCOVER
            success = 1;
        end;
        nepo = 0;
    end;
    if success ==1;
        cvg_n =  (TIME(i,1)* 3600 * 24 * 7 + TIME(i,2)) - (TIME(1,1)* 3600 * 24 * 7 + TIME(1,2));
        rms_n = rms(ENU(i:m,2));
        rms_n0 = rms(ENU(1:i-1,2));
        break;
    end;
end;
if i == m - NCOVER
    rms_n = rms(ENU(1:m,2));
    rms_n0 = rms(ENU(1:m,2));
end;

success = 0;
nepo = 0;
cvg_e = m;
for i=1:m - NCOVER
    if abs(ENU(i,1)) < 0.1
        for j=i:i+NCOVER -1
            if abs(ENU(j,1)) < 0.1
                nepo = nepo + 1;
            end;
        end;
        if nepo == NCOVER
            success = 1;
        end;
        nepo = 0;
    end;
    if success ==1;
        cvg_e =  (TIME(i,1)* 3600 * 24 * 7 + TIME(i,2)) - (TIME(1,1)* 3600 * 24 * 7 + TIME(1,2));
        rms_e = rms(ENU(i:m,1));
        rms_e0 = rms(ENU(1:i-1,1));
        break;
    end;
end;
if i == m -NCOVER
    rms_e = rms(ENU(1:m,1));
    rms_e0 = rms(ENU(1:m,1));
end

success = 0;
nepo = 0;
cvg_u = m;
for i=1:m - NCOVER
    if abs(ENU(i,3)) < 0.1
        for j=i:i + NCOVER - 1
            if abs(ENU(j,3)) < 0.1
                nepo = nepo + 1;
            end;
        end;
        if nepo == NCOVER
            success = 1;
        end;
        nepo = 0;
    end;
    if success ==1;
        cvg_u =  (TIME(i,1)* 3600 * 24 * 7 + TIME(i,2)) - (TIME(1,1)* 3600 * 24 * 7 + TIME(1,2));
        rms_u = rms(ENU(i:m,3));
        rms_u0 = rms(ENU(1:i-1,3));
        break;
    end;
end;
if i == m - NCOVER
    rms_u = rms(ENU(1:m,3));
    rms_u0 = rms(ENU(1:m,3));
end

cvg = max([cvg_e,cvg_n,cvg_u]);

% success = 0;
% nepo = 0;
% cvg = m;
% for i=1:m - NCOVER
%     if abs(ENU(i,1)) < 0.1 && abs(ENU(i,2)) < 0.1 && abs(ENU(i,3)) < 0.1
%         for j=i:i + NCOVER - 1
%             if abs(ENU(j,1)) < 0.1 && abs(ENU(j,2)) < 0.1  && abs(ENU(i,3)) < 0.1
%                 nepo = nepo + 1;
%             end;
%         end;
%         if nepo == NCOVER
%             success = 1;
%         end;
%         nepo = 0;
%     end;
%     if success ==1;
%         cvg =  (TIME(i,1)* 3600 * 24 * 7 + TIME(i,2)) - (TIME(1,1)* 3600 * 24 * 7 + TIME(1,2));
%         break;
%     end;
% end;

output_file=[pathname 'analysis.ana'];
fid=fopen(output_file,'a');
fprintf(fid,'%s',filename(1:4));
fprintf(fid,'  %8.2f',cvg_e/60.0);
fprintf(fid,'  %8.2f',cvg_n/60.0);
fprintf(fid,'  %8.2f',cvg_u/60.0);
% fprintf(fid,'  %8.2f %s',cvg/60.0,'(min)');
fprintf(fid,'  %8.2f',cvg/60.0);
fprintf(fid,'  %10.2f',rms_e0*100);
fprintf(fid,'  %10.3f',rms_n0*100);
fprintf(fid,'  %10.3f',rms_u0*100);
fprintf(fid,'  %10.2f',rms_e*100);
fprintf(fid,'  %10.3f',rms_n*100);
% fprintf(fid,'  %10.3f %s',rms_u*100,'(cm)');
fprintf(fid,'  %10.3f',rms_u*100);
fprintf(fid,'\n');
fclose('all');
% if rms_n <0.2 && rms_e< 0.2 && rms_u < 0.5
%     output_file1=[pathname 'analysis.ana1'];
%     fid1=fopen(output_file1,'a');
%         fprintf(fid1,'%s',filename(1:4));
%         fprintf(fid1,'%10.2f ',cvg_e/60.0);
%         fprintf(fid1,'%10.2f ',cvg_n/60.0);
%         fprintf(fid1,'%10.2f ',cvg_u/60.0);
%         fprintf(fid1,'%10.2f ',cvg/60.0);
%         fprintf(fid1,'%15.2f ',rms_e*100);
%         fprintf(fid1,'%15.3f ',rms_n*100);
%         fprintf(fid1,'%15.3f ',rms_u*100);
%         fprintf(fid1,'\n');
%     fclose('all');
% end;

