function mjd=transT_ymd2mjd(year,month,day)
%SYNTAX:
%     ======================================
%     | mjd=transT_ymd2mjd(year,month,day) |
%     ======================================
%
%     year, month, day or year, day of year to mjd
%
%INPUTS:
%     year: 4-digit year
%    month: month, 1 to 12 or 0 indicates day is day of year
%      day: 1 to 31 or day of year, if imonth is zero
%
%OUTPUT:
%      mjd: integer modified Julian date
%
%originally written by Feng Zhou on 2014/09/22 @ East China Normal Univ.
%
%email: zhouforme@163.com
%
%Dept. of Information Science Technology, East China Normal Univ., China

%%===========================BEGIN PROGRAM====================================%%
doy_of_month=[0,31,59,90,120,151,181,212,243,273,304,334];

% check the input data
if (year<0||month<0||day<0||month>12||day>366)||(month~=0&&day>31)
    error(' *** ERROR (ymd2mjd): incorrect date (year,month,day): %d %d %d', ...
        year,month,day);
end

year1=year;
year1=transT_yy2yyyy(year1);

% doy to month day
if (month==0)
    [im,id]=transT_yrdoy2ymd(year1,day);
else
    im=month;
    id=day;
end
year2=year1;
if (im<=2)
    year2=year2-1;
end
mjd=365*year1-678941+fix(year2/4)-fix(year2/100)+fix(year2/400)+id;
if (im~=0)
    mjd=mjd+doy_of_month(im);
end

return;
%%===========================END PROGRAM======================================%%