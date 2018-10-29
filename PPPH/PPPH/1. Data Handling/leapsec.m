function [leap_sec] = leapsec(mjd)

if mjd<51544
    error('The given date have to be greater than December 1st,2000');
end

if     mjd>=51544 && mjd<53736 %01.01.2000-01.01.2006
    leap_sec = 32;
elseif mjd>=53736 && mjd<54832 %01.01.2006-01.01.2009
    leap_sec = 33;
elseif mjd>=54832 && mjd<56109 %01.01.2009-01.07.2012
    leap_sec = 34;
elseif mjd>=56109 && mjd<57204 %01.07.2012-01.07.2015
    leap_sec = 35;
elseif mjd>=57204 && mjd<57754 %01.07.2015-01.01.2017
    leap_sec = 36;
else
    leap_sec = 37;
end
    
end

