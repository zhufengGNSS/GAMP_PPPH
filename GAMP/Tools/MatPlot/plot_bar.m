function plot_bar(filepath,ana_filename)

%clc;
%clear;
%close all;

%filepath = '\\Mac\Home\Desktop\clkl2_g1\20170101\';
%ana_filename = 'analysis.ana1';

fullname = [filepath  ana_filename];
font_size = 10;

A = importdata(fullname,' ');
B = A.textdata;
C = A.data;
[m,n]=size(C);

%  ’¡≤ ±º‰
figure(1);
bar(C(:,4));
ylabel('Convergence time(min)','FontSize',15);
set(gca,'XTick',1:1:m);
set(gca,'XTickLabel',B);
h = gca;
th = rotateticklabel(h,90,font_size);
xlim([0 m+1]);
saveas(gcf,[filepath  'convergence.jpg'],'jpg');
saveas(gcf,[filepath  'convergence.fig'],'fig');
delete(gcf);

figure(2);
bar(C(:,5));
ylabel('North RMS(cm)','FontSize',15);
set(gca,'XTick',1:1:m);
set(gca,'XTickLabel',B);
h = gca;
th = rotateticklabel(h,90,font_size);
xlim([0 m+1]);
saveas(gcf,[filepath  'RMS_N.jpg'],'jpg');
saveas(gcf,[filepath  'RMS_N.fig'],'fig');
delete(gcf);

figure(3);
bar(C(:,6));
ylabel('East RMS(cm)','FontSize',15);
set(gca,'XTick',1:1:m);
set(gca,'XTickLabel',B);
h = gca;
th = rotateticklabel(h,90,font_size);
xlim([0 m+1]);
saveas(gcf,[filepath  'RMS_E.jpg'],'jpg');
saveas(gcf,[filepath  'RMS_E.fig'],'fig');
delete(gcf);

figure(4);
bar(C(:,7));
ylabel('Up RMS(cm)','FontSize',15);
set(gca,'XTick',1:1:m);
set(gca,'XTickLabel',B);
h = gca;
th = rotateticklabel(h,90,font_size);
xlim([0 m+1]);
saveas(gcf,[filepath  'RMS_U.jpg'],'jpg');
saveas(gcf,[filepath  'RMS_U.fig'],'fig');
delete(gcf);
