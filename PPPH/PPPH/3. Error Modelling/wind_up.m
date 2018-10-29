function [wup] = wind_up(rec,sat,sun,prev)

esun = sun - sat;
esun = esun./norm(esun);
ez   = -1.*sat; ez = ez./norm(ez);
ey   = cross(ez,esun); ey = ey./norm(ey);
ex   = cross(ey,ez);   ex = ex./norm(ex);
xs   = ex; ys = ey;

[elip] = xyz2plh(rec,0);
phi = elip(1);
lam = elip(2);

xr   = [(-sin(phi)*cos(lam)) (-sin(phi)*sin(lam)) cos(phi)];
yr   = [  sin(lam) -cos(lam) 0];

k = rec - sat; k = k./norm(k);

Ds = xs - k.*(dot(k,xs)) - cross(k,ys);
Dr = xr - k.*(dot(k,xr)) + cross(k,yr);
wup= acos(dot(Ds,Dr)/norm(Ds)/norm(Dr));
if dot(k,(cross(Ds,Dr)))<0, wup = -wup; end
wup = (2*pi*floor(((prev - wup)/2/pi)+0.5)) + wup;
end

