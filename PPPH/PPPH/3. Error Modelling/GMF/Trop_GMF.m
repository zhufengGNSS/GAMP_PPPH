function [Trop,Mwet,Mn,Me,ZHD] = Trop_GMF(rec,sat,dmjd,p)

[ellp] = xyz2plh(rec,0);
dlat = ellp(1);
dlon = ellp(2);
hell = ellp(3);

[Az,Elv] = local(rec,sat,0);

f = 0.0022768;
k = 1 - (0.00266*cos(2*dlat)) - (0.28*10^-6*hell);
ZHD = f*(p/k);

[gmfh,gmfw] = gmf_f_hu(dmjd,dlat,dlon,hell,(pi/2 - Elv));

Trop = gmfh*ZHD;
Mwet = gmfw;
Mg = 1/((tan(Elv)*sin(Elv))+0.0032);
Mn = Mg*cos(Az);
Me = Mg*sin(Az);
end

