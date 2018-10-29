function plot_neu(pathname,filename,mode,ylim1)
mark = 5;
width = 1.0;
staname=strrep(filename,'_','+');
fullname = [pathname  filename];
A=load(fullname);
%gpst = mod(A(:,2),86400*2)/3600;0;
diff_week = A(:,7)-A(1,7);  %与第一个历元的gpsweek之差
gpst = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
X=A(:,12);
Y=A(:,13);
Z=A(:,14);
figure;
plot(gpst,X,'b-','markersize',mark,'LineWidth',width);
hold on;
plot(gpst,Y,'-','Color',[0,127,0]/255,'markersize',mark,'LineWidth',width);
hold on;
plot(gpst,Z,'r-','markersize',mark,'LineWidth',width);
hold on;
% plot(gpst,Z,'b-','markersize',mark,'LineWidth',width);
% hold on;
% plot(gpst,X,'g-','markersize',mark,'LineWidth',width);
% hold on;
% plot(gpst,Y,'r-','markersize',mark,'LineWidth',width);
% hold on;

title(staname);
%len = length(gpst);
%xlim([0,gpst(len)*1.005]);
if ylim1 ~=0
    ylim([-1*ylim1,ylim1]);
    set(gca,'ytick',-1*ylim1:ylim1/10.0:ylim1);
end;
xlabel('Epoch [hour]');
ylabel('Error [m]');
legend('E','N','U');
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
% pathname_fig = sprintf('%s%s',pathname,'plot\fig\');
% if ~exist(pathname_fig) 
%     mkdir(pathname_fig)
% end 
% saveas(gcf,[pathname_fig staname '.fig'],'fig');
delete(gcf);