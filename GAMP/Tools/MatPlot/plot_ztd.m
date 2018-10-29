function plot_ztd(pathname,filename,mode,ylim1)
staname=strrep(filename,'_','+');
fullname = [pathname  filename];
A=load(fullname);
diff_week = A(:,7)-A(1,7);  %与第一个历元的gpsweek之差
gpst = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
X=A(:,11);
% intevel = 1;
% s_epo = intevel;
% e_epo = length(X)/intevel;
figure;
% plot(s_epo:intevel:e_epo,X,'-');
plot(gpst,X,'-','LineWidth',2);
hold on;

%s_epo = A(1,4)+ A(1,5)/60;
%e_epo = A(length(X),4)+ A(length(X),5)/60;
%intevel = 1/12.0;

%set(gca,'Xtick',[round(s_epo):1:round(e_epo)]);
%set(gca,'Xtick',[round(s_epo):1:round(e_epo)]);
% set(gcf,'visible','off');
title(staname);
box on;
if ylim1 ~=0
    ylim([-1*ylim1,ylim1]);
    set(gca,'ytick',-1*ylim1:ylim1/10.0:ylim1);
end;
xlabel('Epoch [hour]');
ylabel('ZTD [m]');
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