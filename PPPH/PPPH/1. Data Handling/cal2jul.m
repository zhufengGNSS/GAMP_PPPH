function [jd,mjd] = cal2jul(year,mon,day,sec)

narginchk(4,4)

sec = sec/3600;

if ~isscalar(year) || ~isscalar(mon) || ~isscalar(day)
    error('Year, Month and Day should be scalar.')
end

if mon<1 || mon>12
    error('Month should be between 1 and 12.')
end


if day<1 || day>31
    error('Day should be between 1 and 31.')
end

if mon<=2
    m = mon + 12;
    y = year - 1;
else
    m = mon;
    y = year;
end

jd = floor(365.25*y) + floor(30.6001*(m+1)) + day + (sec/24) + 1720981.5;
mjd= jd - 2400000.5;

nargoutchk(1,2)
end

