function plot_wlamb(pathname,filename,mode)
staname=strrep(filename,'_','+');
fullname = [pathname  filename];
A=load(fullname);
A(find(A == 99999.0)) = NaN;
%gpst = mod(A(:,2),86400*2)/3600;0;
diff_week = A(:,7)-A(1,7);  %与第一个历元的gpsweek之差
gpst = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
X=A(:,9:size(A,2));
figure;
plot(gpst,X);
hold on;
% ylim([-20 20]);
title(staname);
xlabel('Epoch [hour]');
ylabel('MW_AMB [cycles]');
box on;
grid on;
hold off;


%index_dir=findstr(pathname,'\');
%path_temp=pathname(1:index_dir(end)-1);
%pathname_jpg = sprintf('%s%s%s',path_temp,'_plot\jpg\');
pathname_jpg = sprintf('%s%s',pathname,'plot\jpg\');
if ~exist(pathname_jpg) 
    mkdir(pathname_jpg)
end 
saveas(gcf,[pathname_jpg staname '.jpg'],'jpg');
%pathname_fig = sprintf('%s%s%s',path_temp,'_plot\fig\');
pathname_fig = sprintf('%s%s',pathname,'plot\fig\');
if ~exist(pathname_fig) 
    mkdir(pathname_fig)
end 
saveas(gcf,[pathname_fig staname '.fig'],'fig');
delete(gcf);