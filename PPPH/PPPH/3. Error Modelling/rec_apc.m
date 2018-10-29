function [rapc] = rec_apc(s_xyz,r_xyz,r_apc,k)

l   = r_xyz - s_xyz;
los = l./norm(l);

[elip] = xyz2plh(r_xyz,0);
lat = elip(1);
lon = elip(2);

ori = [(-sin(lon)) (-cos(lon)*sin(lat)) (cos(lon)*cos(lat));...
        (cos(lon)) (-sin(lon)*sin(lat)) (sin(lon)*cos(lat));...
               (0) (cos(lat)) (sin(lat))];
           
f = r_apc(:,:,k);

c = [f(2);
     f(1);
     f(3)];

p = ori*c;
if size(los,1)~=size(p,1)
    los = los';
end

rapc = dot(p,los);

end

