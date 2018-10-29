path = 'F:\output\ztd_result20161021\';
path1 = sprintf('%s%s%s',path,'20161021fix_5','\');
path2 = sprintf('%s%s%s',path,'20161021fix_6','\');
path3 = sprintf('%s%s%s',path,'20161021fix_7','\');
path4 = sprintf('%s%s%s',path,'20161021fix_8','\');
path5 = sprintf('%s%s%s',path,'igs_ztd','\');

filenames1 = dir([path1 '*.ztd']);
filenames2 = dir([path2 '*.ztd']);
filenames3 = dir([path3 '*.ztd']);
filenames4 = dir([path4 '*.ztd']);
filenames5 = dir([path5 '*.ztd']);

index1 = -1;
index2 = -1;
index3 = -1;
index4 = -1;
index5 = -1;

for i = 1 : length(filenames1)
    filename1 = filenames1(i,1).('name');
    staname1 = filename1(1:4);
    fullname1 = [path1  filename1];
    index1 = i;
    for j = 1 : length(filenames2)
        filename2 = filenames2(j,1).('name');
        staname2 = filename2(1:4);
        fullname2 = [path2  filename2];
        if staname2 ==staname1
            index2 =j;
            break;
        end;
    end;
    for j = 1 : length(filenames3)
        filename3 = filenames3(j,1).('name');
        staname3 = filename3(1:4);                
        fullname3 = [path3  filename3];
        if staname3 ==staname1
            index3 =j;
            break;
        end;
    end;    
    for j = 1 : length(filenames4)
        filename4 = filenames4(j,1).('name');
        staname4 = filename4(1:4);
        fullname4 = [path4  filename4];
        if staname4 ==staname1
            index4 =j;
            break;
        end;
    end;
    for j = 1 : length(filenames5)
        filename5 = filenames5(j,1).('name');
        staname5 = filename5(1:4);
        fullname5 = [path5  filename5];
        if staname5 ==staname1
            index5 =j;
            break;
        end;
    end;
    if index2==-1||index3==-1 ||index4==-1 ||index5==-1
        continue;
    end;
    A1=load(fullname1);
    A2=load(fullname2);
    A3=load(fullname3);
    A4=load(fullname4);
    A5=load(fullname5);
    
    figure(i);
    diff_week = A1(:,7)-A1(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A1(:,8)/3600;% - A1(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
    ZTD1=A1(:,9);
    %plot(gpst,ZTD1,'g.');
    hold on;

    diff_week = A2(:,7)-A2(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A2(:,8)/3600;%  - A2(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
    ZTD2=A2(:,9);
    plot(gpst,ZTD2,'r.');
    hold on;
    
    diff_week = A3(:,7)-A3(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A3(:,8)/3600;%  - A3(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
    ZTD3=A3(:,9);
    plot(gpst,ZTD3,'y.');
    hold on;
    
    diff_week = A4(:,7)-A4(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A4(:,8)/3600;%  - A4(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
    ZTD4=A4(:,9);
    %plot(gpst,ZTD4,'b.');
    hold on;
    
    diff_week = A5(:,7)-A5(1,7);  %与第一个历元的gpsweek之差
    gpst = diff_week *7 *24 + A5(:,8)/3600;%  - A5(1,8)/3600;  %当前历元在周内的小时数 - 第一个历元的小时数
    ZTD5=A5(:,9);
    plot(gpst,ZTD5,'k.');
    hold on;
    
    grid on;
    box on;
    title(staname1);
    xlabel('Epoch(s)');
    ylabel('ZTD(mm)');
    %legend('1E+5','1E+6','1E+7','1E+8','igs');
    legend('1E+6','1E+7','igs');

    hold off;
    saveas(figure(i),[path1 staname1 '.jpg'],'jpg');
    delete(gcf);
end;
