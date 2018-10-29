
function plot_panda()
global newpath;
[filenames,pathname] = uigetfile('*','open nsat files','MultiSelect','on',newpath);
n= max(size(filenames));
if ~iscell(filenames)
    plot_neu2(pathname , filenames);
else
    for i=1:n
        filename=filenames{i};
        plot_neu2(pathname , filename);
    end
end;
newpath = pathname;


function plot_neu2(pathname , filename)
mark = 5;
width = 2.5;
ylim1 = 1;
fullname = [pathname filename];
A = load(fullname,' ');
gpst =  A(:,1)/3600 - A(1,1)/3600;
X=A(:,2);
Y=A(:,3);
Z=A(:,4);
figure;
plot(gpst,Z,'k-','markersize',mark,'LineWidth',width);
hold on;
plot(gpst,Y,'r-','markersize',mark,'LineWidth',width);
hold on;
plot(gpst,X,'g-','markersize',mark,'LineWidth',width);
hold on;
ylim([-1*ylim1,ylim1]);
set(gca,'ytick',-1*ylim1:ylim1/10.0:ylim1);
grid on;
title(filename(13:16));
legend('U','N','E');
saveas(gcf,[pathname filename(13:16) '.jpg'],'jpg');
delete(gcf);