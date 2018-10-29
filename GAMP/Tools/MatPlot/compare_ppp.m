ylim1 = 1.0;
staname = 'CUT0';
path = sprintf('%s%s%s','F:\output\data2\',staname,'\');
filename1 = sprintf('%s%s',path,staname,'_g.ppp');
filename2 = sprintf('%s%s',path,staname,'_grec.ppp');
A=load(filename1);
B=load(filename2);
%gpst = mod(A(:,2),86400*2)/3600;0;
diff_week = A(:,7)-A(1,7);  %与第一个历元的gpsweek之差
gpst1 = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
X1=A(:,9);
Y1=A(:,10);
Z1=A(:,11);
diff_week2 = B(:,7)-B(1,7);  %与第一个历元的gpsweek之差
gpst2 = diff_week2 *7 *24 + B(:,8)/3600 - B(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
X2=B(:,9);
Y2=B(:,10);
Z2=B(:,11);

figure(1);
%subplot(3,1,1);
plot(gpst1,X1,'g.',gpst2,X2,'k.');
grid on;
title([staname '(N)']);
ylim([ylim1*-1,ylim1]);
xlabel('Epoch(hour)');
ylabel('Error(m)');
legend('G','G/R/E/C');
saveas(gcf,[path staname '+N.fig'],'fig');
saveas(gcf,[path staname '+N.jpg'],'jpg');

figure(2);
%subplot(3,1,2);
plot(gpst1,Y1,'g.',gpst2,Y2,'k.');
grid on;
title([staname '(E)']);
ylim([ylim1*-1,ylim1]);
xlabel('Epoch(hour)');
ylabel('Error(m)');
legend('G','G/R/E/C');
saveas(gcf,[path staname '+E.fig'],'fig');
saveas(gcf,[path staname '+E.jpg'],'jpg');

figure(3);
%subplot(3,1,3);
plot(gpst1,Z1,'g.',gpst2,Z2,'k.');
grid on;
title([staname '(U)']);
ylim([ylim1*-1,ylim1]);
xlabel('Epoch(hour)');
ylabel('Error(m)');
legend('G','G/R/E/C');
saveas(gcf,[path staname '+U.fig'],'fig');
saveas(gcf,[path staname '+U.jpg'],'jpg');