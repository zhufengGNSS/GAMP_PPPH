function [data] = cs_detect(data,options)

[arc]  = arc_dtr(data.obs);

sn = size(data.obs.st,2);

c  = 299792458; % m/s
[freq,wavl] = frequencies;

dt  = data.inf.time.int;
sig0 = sqrt(2*(0.0027^2 + 0.0017^2));%meter
dl  = (0.4)*(dt/3600); %meter/hour

for k=1:sn
    f1 = freq(k,1); f2 = freq(k,2);
    lamwl = c/(f1-f2);
    
    lwl = (data.obs.l1(:,k).*f1 - data.obs.l2(:,k).*f2)./(f1-f2);
    pnl = (data.obs.p1(:,k).*f1 + data.obs.p2(:,k).*f2)./(f1+f2);
    nwl = (lwl - pnl)./(lamwl);
    
    gfl = data.obs.l1(:,k) - data.obs.l2(:,k); %meter
    
    ark = arc{k};
    for t=1:size(ark,1)
        st = ark(t,1);
        fn = ark(t,2);
        
        for i=(ark(t,1)+1):ark(t,2)
            
            dmwc = 0;
            dgfc = 0;
            
            if options.CSMw == 1
                if (i-2<st)
                    mmw = mean(nwl(st:i-1,1));
                    smw =  std(nwl(st:i  ,1));
                elseif (i-30<st)
                    mmw = mean(nwl(st:i-1,1));
                    smw =  std(nwl(st:i-1,1));
                else
                    mmw = mean(nwl(i-30:i-1,1));
                    smw =  std(nwl(i-30:i-1,1));
                end
                dmw = mmw - nwl(i,1);
                if abs(dmw)>(5*smw)
                    dmwc = 1; 
                end
            end
            
            if options.CSGf == 1
                dgf = gfl(i-1,1) - gfl(i,1);
                elv = data.obs.elv(i,k);
                me = 1 + (10*exp(-elv/10));
                smg = sig0*me;
                
                if abs(dgf)>((4*smg)+dl)
                    dgfc = 1; 
                end
            end
            
            if (dmwc==1 && std(nwl(st:fn,1))>0.6) || dgfc==1
                one  = nwl(i-1,1) - nwl(i,1);
                two  = gfl(i-1,1) - gfl(i,1);
                A = [1 -1;wavl(k,1) -wavl(k,1)];
                L = [one;two];
                Dn = pinv(A)*L;
                Dn1 = round(Dn(1));
                Dn2 = round(Dn(2));
                
                if (Dn1~=0) && (Dn2~=0)
                    data.obs.l1(i:fn,k) = data.obs.l1(i:fn,k) + Dn1.*wavl(k,1);
                    data.obs.l2(i:fn,k) = data.obs.l2(i:fn,k) + Dn2.*wavl(k,2);
                end
                
                st = i;
            end
        end
    end
end
end