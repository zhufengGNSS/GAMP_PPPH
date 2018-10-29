path = 'F:\output\ztd_result20161021\20161021fix_7\com\';
filenames = dir([path '*.comztd']);
for i = 1 : length(filenames)    
    filename = filenames(i,1).('name');
    staname = filename(1:4);
    fullname = [path  filename];
    A=load(fullname);
    diff_ztd = A(:,10)-A(:,9);
    rms = std(diff_ztd);
    diff_week = A(:,7)-A(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A(:,8)/3600 - A(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数    
    figure(i);
    plot(gpst,A(:,9),'r.',gpst,A(:,10),'b.');
    grid on;
    box on;
    titlename = sprintf('%s %s%.2f%s',staname,'(rms:',rms,'mm)');
    title(titlename);
    xlabel('Epoch(s)');
    ylabel('ZTD(mm)');
    legend('rtppp','igs');
    saveas(figure(i),[path staname '.jpg'],'jpg');
    saveas(figure(i),[path staname '.fig'],'fig');
    delete(gcf);
end;