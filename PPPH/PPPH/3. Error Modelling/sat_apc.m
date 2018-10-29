function [sapc] = sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,opt,sno)

l   = r_xyz - s_xyz;
los = l./norm(l);
% satellite -fixed coordinate system
k = (-1).*(s_xyz./(norm(s_xyz)));
rs  = sun_xyz - s_xyz;
e   = rs./(norm(rs));
j   = cross(k,e);
i   = cross(j,k);
sf  = [i; j; k];

de1 = (s_apc(:,:,1))';
de2 = (s_apc(:,:,2))';

if (any(isnan(de1)) || any(isnan(de2))) && (sno>58 && sno<89) % GALILEO CONVENTIONAL
    de1 = [0.2;0;0.6];
    de2 = [0.2;0;0.6];
elseif (any(isnan(de1)) || any(isnan(de2))) && (sno>88 && sno<93) % BEIDOU CONVENTIONAL
    de1 = [0.6;0;1.1];
    de2 = [0.6;0;1.1];
end
    
if opt == 1
    rk = sf\de1;
elseif opt == 2
    rk = sf\de2;
end

sapc = dot(rk,los);
end