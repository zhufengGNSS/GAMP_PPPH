function [data] = outlier(data)

c  = 299792458; % m/s
[freq,~] = frequencies;

[arc]  = arc_dtr(data.obs);


sn = size(data.obs.st,2);


for k=1:sn
    
    f1 = freq(k,1); f2 = freq(k,2);
    
    lwl = (data.obs.l1(:,k).*f1 - data.obs.l2(:,k).*f2)./(f1-f2);
    pnl = (data.obs.p1(:,k).*f1 + data.obs.p2(:,k).*f2)./(f1+f2);
    lamwl = c/(f1 - f2);
    nwl = (lwl - pnl)./(lamwl);
    
    ark = arc{k};
    
    for n=1:size(ark,1)
        st = ark(n,1);
        fn = ark(n,2);
        while 1
            t = find(data.obs.st(:,k)==1);
            t = t(t(:)>=st & t(:)<=fn);
            L = nwl(t);
            ran = size(t,1);
            A = [(t).^2 t ones(ran,1)];
            
            X = A\L;
            V = L - A*X;
            
            rmse = sqrt(sum(V.^2)/ran);
            det = find(abs(V)>(4*0.6));
            if rmse>0.6 && any(det)
                
                del = t(abs(V)==max(abs(V)));
                data.obs.st(del,k) = 0;
            else
                break
            end
        end
    end
end
end

