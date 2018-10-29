station = 'NCKU';
mark = 5;
width = 2.5;
fpath = 'C:\gnssdata\lixin\';
A = load( 'C:\gnssdata\lixin\lixin_par_result\par_183\enu_2016183_ncku');
B = load( ['C:\gnssdata\lixin\mgex\183\res_c_10\20160701\' station '_20160701_sp3_C.ppp']);
C = load( ['C:\gnssdata\lixin\mgex\183\res_c_1\20160701\' station '_20160701_sp3_C.ppp']);
D = load( ['C:\gnssdata\lixin\mgex\183\res_c_0.5\20160701\' station '_20160701_sp3_C.ppp']);
%E = load( [fpath station '_20161120_sp3_GC.ppp']);
%F = load( [fpath station '_20161120_sp3_GRC.ppp']);

%figure;    

n = 1;
for i=9:11    
    figure;
    diff_week = B(:,7)-B(1,7);
    gpst = diff_week *7 *24 + B(:,8)/3600 - B(1,8)/3600;
    plot(gpst,B(:,i),'g-','markersize',mark,'LineWidth',width);
    hold on;
    diff_week = C(:,7)-C(1,7);
    gpst = diff_week *7 *24 + C(:,8)/3600 - C(1,8)/3600;
    plot(gpst,C(:,i),'c-','markersize',mark,'LineWidth',width);
    hold on;
    diff_week = D(:,7)-D(1,7);
    gpst = diff_week *7 *24 + D(:,8)/3600 - D(1,8)/3600;
    plot(gpst,D(:,i),'r-','markersize',mark,'LineWidth',width);
    hold on;
    %diff_week = E(:,7)-E(1,7);
    %gpst = diff_week *7 *24 + E(:,8)/3600 - E(1,8)/3600;
    %plot(gpst,E(:,i),'m-','markersize',mark,'LineWidth',width);
    %hold on;
    %diff_week = F(:,7)-F(1,7);
    %gpst = diff_week *7 *24 + F(:,8)/3600 - F(1,8)/3600;
    %plot(gpst,F(:,i),'y-','markersize',mark,'LineWidth',width);
    gpst = A(:,1)/3600;
    plot(gpst,A(:,i-7),'k-','markersize',mark,'LineWidth',width);
    hold on;
    if i>10
        lim =1;
    else
        lim =0.5;
    end;
    ylim([lim * (-1),lim]);
    if i ==9
        de1 = '+N';
        re1 = '(North)';
    elseif i ==10
        de1 = '+E';
        re1 = '(East)';
    elseif i ==11
        de1 = '+U';
        re1 = '(Up)';
    end;
    xlabel('Epoch(hour)');
    ylabel('Error(m)');
    %legend('G','R','C','GR','GC','GRC');
    legend('10.0','1.0','0.5','panda');
    title( [station re1]);
    grid on;
    saveas(gcf,[fpath station de1 '.jpg'],'jpg');
end;