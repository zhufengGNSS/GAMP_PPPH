station = 'CUT1';
mark = 5;
width = 2.5;
fpath = 'C:\glonass\cut0\20161120\';

A = load( [fpath station '_20161120_sp3_G.ppp']);
B = load( [fpath station '_20161120_sp3_R.ppp']);
C = load( [fpath station '_20161120_sp3_C.ppp']);
D = load( [fpath station '_20161120_sp3_GR.ppp']);
E = load( [fpath station '_20161120_sp3_GC.ppp']);
F = load( [fpath station '_20161120_sp3_GRC.ppp']);

%figure;    

n = 1;
for i=9:11    
    %subplot(3,1,n);
    %n = n + 1;
    figure;
    diff_week = A(:,7)-A(1,7);
    gpst = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;
    plot(gpst,A(:,i),'r-','markersize',mark,'LineWidth',width);
    hold on;
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
    plot(gpst,D(:,i),'y-','markersize',mark,'LineWidth',width);
    hold on;
    diff_week = E(:,7)-E(1,7);
    gpst = diff_week *7 *24 + E(:,8)/3600 - E(1,8)/3600;
    plot(gpst,E(:,i),'m-','markersize',mark,'LineWidth',width);
    hold on;
    diff_week = F(:,7)-F(1,7);
    gpst = diff_week *7 *24 + F(:,8)/3600 - F(1,8)/3600;
    plot(gpst,F(:,i),'k-','markersize',mark,'LineWidth',width);
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
    legend('G','R','C','GR','GC','GRC');
    title( [station re1]);
    grid on;
    saveas(gcf,[fpath station de1 '.jpg'],'jpg');
end;