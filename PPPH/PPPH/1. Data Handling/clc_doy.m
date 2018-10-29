function [doy] = clc_doy(year,mon,day)

if rem(year,400)==0
    sit = 1;
elseif rem(year,300)==0
    sit = 0;
elseif rem(year,200)==0
    sit = 0;
elseif rem(year,100)==0
    sit = 0;
elseif rem(year,4)==0
    sit = 1;
else
    sit = 0;
end

switch sit
    case 0
        dom = [31 28 31 30 31 30 31 31 30 31 30 31];
    case 1
        dom = [31 29 31 30 31 30 31 31 30 31 30 31];
end

if mon-1<1
    doy = day;
else
    doy = sum(dom(1:(mon-1))) + day;
end

end

