function [azim,elev] = local(rec,sat,dopt)

narginchk(3,3)

if numel(rec)~=3 && numel(sat)~=3
    error('Receiver and satellite position vectors must include X,Y,Z')
end

if size(rec,1)~=size(sat,1)
    if size(rec,1)==3
        rec = rec';
    elseif size(sat,1)==3
        sat = sat';
    end
end

los = sat - rec;


p   = los./norm(los);


[ellp] = xyz2plh(rec,0);
lat = ellp(1);
lon = ellp(2);


e = [-sin(lon)           cos(lon)                 0];
n = [-cos(lon)*sin(lat) -sin(lon)*sin(lat) cos(lat)];
u = [ cos(lon)*cos(lat)  sin(lon)*cos(lat) sin(lat)];


elev = asin(dot(p,u));
azim = atan2(dot(p,e),dot(p,n));
azim = mod(azim,2*pi);

if dopt == 1
    t = (180/pi);
    elev = elev.*t;
    azim = azim.*t;
end

nargoutchk(1,2)
end

