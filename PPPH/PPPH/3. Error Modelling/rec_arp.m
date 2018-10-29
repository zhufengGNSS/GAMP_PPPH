function [rarp] = rec_arp(s_xyz,r_xyz,arp)

l   = r_xyz -s_xyz;
los = l./norm(l);

[elip] = xyz2plh(r_xyz,0);
lat = elip(1);
lon = elip(2);

ori = [(-sin(lon)) (-cos(lon)*sin(lat)) (cos(lon)*cos(lat));...
        (cos(lon)) (-sin(lon)*sin(lat)) (sin(lon)*cos(lat));...
               (0)           (cos(lat))          (sin(lat))];
           
enu = [arp(2);arp(3);arp(1)];
p = ori*enu;
if size(p,1)~=size(los,1)
    los = los';
end
rarp = dot(p,los);
end

