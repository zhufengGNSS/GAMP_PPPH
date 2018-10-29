function [model] = nmodel(data,options)


c = 299792458;

[~,wavl] = frequencies;

mopt = zeros(1,11);
mopt( 1) = options.SatClk;
mopt( 2) = options.SatAPC;
mopt( 3) = options.RecAPC;
mopt( 4) = options.RecARP;
mopt( 5) = options.RelClk;
mopt( 6) = options.SatWind;
mopt( 7) = options.AtmTrop;
mopt( 9) = options.RelPath;
mopt(11) = options.Solid;

satno  = dtr_satno(data.obs);
n_sats = sum(satno);
n_model= n_sats*4;
model  = zeros(n_model,31);

if strcmp(data.inf.time.system,'GPS') || isempty(data.inf.time.system)
    dt = 51.184;
else
    dt  = 32.184 + leap;
end

r_xyz = data.inf.rec.pos;
if size(r_xyz,1)~=1
    r_xyz = r_xyz';
end

[ellp] = xyz2plh(r_xyz,0);
dlat = ellp(1);
dlon = ellp(2);
hell = ellp(3);
[~,dmjd] = cal2jul(data.inf.time.first(1),data.inf.time.first(2),data.inf.time.first(3),0);

[p] = ngpt2(dmjd,dlat,dlon,hell,1,1);

arp = data.inf.ant.hen;
if size(arp,1)~=1
    arp = arp';
end

r_apc = data.atx.rcv.neu;


t = 0;

for i=1:size(data.obs.st,1)
    
    year = data.inf.time.first(1); 
    doy = data.inf.time.doy; 
    secod = data.obs.ep(i,1);
    
    [~,mjd_tt] = cal2jul(data.inf.time.first(1),data.inf.time.first(2),data.inf.time.first(3),...
                 (secod+dt));
             
    [~,mjd] = cal2jul(data.inf.time.first(1),data.inf.time.first(2),data.inf.time.first(3),...
                 (secod)); 
             
    sun_xyz = sun(mjd_tt); 
    mon_xyz = moon(mjd_tt);
    sats = find(data.obs.st(i,:)==1);
    for k = sats
        
        s_xyz = data.psat(i,1:3,k);
        v_xyz = data.psat(i,4:6,k);
        s_apc = data.atx.sat.neu(k,:,:);
        for u=1:4
            t = t + 1;
            
            model(t,1) = year;
            
            model(t,2) = doy;
            
            model(t,3) = secod;
            
            model(t,4) = k;
            
            model(t,5) = data.tofs(i,k);
            
            switch u
                case 1 
                    model(t,6 ) = data.obs.p1(i,k);
                    
                    sdt = -(sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,1,k));
                    if ~isnan(sdt)
                        model(t,16) = sdt;
                    end
                    
                    if k>32 && k<59
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,3);
                    else
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,1);
                    end
                case 2
                    model(t,6) = data.obs.p2(i,k);
                    
                    sdt = -(sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,2,k));
                    if ~isnan(sdt)
                        model(t,16) = sdt;
                    end
                    
                    if k>32 && k<59
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,4);
                    else
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,2);
                    end
                case 3
                    
                    model(t,6) = data.obs.l1(i,k);
                    
                    sdt = -(sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,1,k));
                    if ~isnan(sdt)
                        model(t,16) = sdt;
                    end
                    if k>32 && k<59
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,3);
                    else
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,1);
                    end
                    
                    if i==1
                        prev=0;
                    else
                        prev = model(t-4,22);
                    end
                    model(t,20) = ((wind_up(r_xyz,s_xyz,sun_xyz,prev))/2/pi)*wavl(k,1);
                case 4
                    model(t,6) = data.obs.l2(i,k);
                    sdt = -(sat_apc(s_xyz,r_xyz,sun_xyz,s_apc,2,k));
                    if ~isnan(sdt)
                        model(t,16) = sdt;
                    end
                    
                    if k>32 && k<59
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,4);
                    else
                        model(t,17) = rec_apc(s_xyz,r_xyz,r_apc,2);
                    end
                    
                    if i==1
                        prev=0;
                    else
                        prev = model(t-4,22);
                    end
                    model(t,20) = ((wind_up(r_xyz,s_xyz,sun_xyz,prev))/2/pi)*wavl(k,1);
            end
            
            model(t,8:10)  = s_xyz;
            
            model(t,11:13) = v_xyz;
            
            model(t,14) = norm(s_xyz - r_xyz);
            
            model(t,15) = -(data.psat(i,7,k)*c);
            
            model(t,18) = rec_arp(s_xyz,r_xyz,arp);
            
            model(t,19) = -(rel_clk(s_xyz,v_xyz));
            
            [Trop,Mwet,Mn,Me,ZHD] = Trop_GMF(r_xyz,s_xyz,mjd,p);
            model(t,21) = Trop;
            model(t,28) = Mwet;
            model(t,29) = Mn;
            model(t,30) = Me;
            model(t,31) = ZHD;
            
            model(t,23) = rpath(r_xyz,s_xyz);
            
            model(t,25) = solid(r_xyz,s_xyz,sun_xyz,mon_xyz);
            
            model(t,26) = data.obs.elv(i,k);
            
            model(t,27) = data.obs.azm(i,k);
            
            
            full = model(t,15:25);
            model(t,7) = sum(full(mopt==1));
        end
    end
end
end

