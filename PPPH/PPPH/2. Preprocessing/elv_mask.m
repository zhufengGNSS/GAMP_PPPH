function [data] = elv_mask(data,options)

en = size(data.obs.st,1);
sn = size(data.obs.st,2);

r_xyz = data.inf.rec.pos;
elv = NaN(en,sn);
azm = NaN(en,sn);
for i=1:en
    for k=1:sn
        if data.obs.st(i,k) == 1
            
            s_xyz = data.psat(i,1:3,k);
            
            [az,elev] = local(r_xyz,s_xyz,1);
            elv(i,k) = elev;
            azm(i,k) = az;
            if elev<options.elvangle
                data.obs.st(i,k) = 0;
            end
        end
    end
end

data.obs.elv = elv;
data.obs.azm = azm;
end

