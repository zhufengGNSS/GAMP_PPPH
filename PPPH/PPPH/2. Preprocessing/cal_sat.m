function [data] = cal_sat(data)

c = 299792458;        % m/s velocity of light
we=7.2921151467e-5;   % WGS84 value for Earth angular velocity (rad/s)

clkint = data.inf.time.clkint;
sp3int = data.inf.time.sp3int;
rec    = data.inf.rec.pos;
if size(rec,1) ~= 1
    rec = rec';
end

en = size(data.obs.st,1);
sn = size(data.obs.st,2);
psat = NaN(en,7,sn);
tofs = NaN(en,sn);


for i=1:en
    for k=1:sn
        if data.obs.st(i,k) == 1
            
            if any(data.obs.p1(i,k))
                tof = data.obs.p1(i,k)/c;
            elseif any(data.obs.p2(i,k))
                tof = data.obs.p2(i,k)/c;
            end
            
            nep = data.obs.ep(i,1) - tof;
            d_clk = data.clk(:,k);
            dt  = entrp(nep,clkint,d_clk);
            
            if isnan(dt)
                data.obs.st(i,k) = 0;
                continue
            end
            nep = nep - dt;
            
            tofs(i,k) = tof + dt;
            
            d_x = data.sat.sp3(:,1,k);
            d_y = data.sat.sp3(:,2,k);
            d_z = data.sat.sp3(:,3,k);
            if data.opt.entrp == 1
                [X,VX] = entrp_orbt(nep,sp3int,d_x);
                [Y,VY] = entrp_orbt(nep,sp3int,d_y);
                [Z,VZ] = entrp_orbt(nep,sp3int,d_z);
            else
                [X,VX] = entrp(nep,sp3int,d_x);
                [Y,VY] = entrp(nep,sp3int,d_y);
                [Z,VZ] = entrp(nep,sp3int,d_z);
            end
            d_clk = data.clk(:,k);
            [dt,~] = entrp(nep,clkint,d_clk);
            
            R  = [X Y Z];
            
            tf     = norm(R - rec)./c;
            er_ang = rad2deg(tf*we);
            pos    = rotation(R,er_ang,3);
            
            if any(isnan(pos))
                data.obs.st(i,k) = 0;
            else
                psat(i,1,k) = pos(1);
                psat(i,2,k) = pos(2);
                psat(i,3,k) = pos(3);
                psat(i,4,k) = VX;
                psat(i,5,k) = VY;
                psat(i,6,k) = VZ;
                psat(i,7,k) = dt;
            end
        end
    end
end

data.psat = psat;
data.tofs = tofs;
end

