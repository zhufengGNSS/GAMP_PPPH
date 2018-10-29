function [data] = smoothing(data)

[freq,~] = frequencies;

[arc] = arc_dtr(data.obs);

sn = size(data.obs.st,2);

for k=1:sn
    f1 = freq(k,1); f2 = freq(k,2);
    
    ark = arc{k};
    for t=1:size(ark,1)
        st = ark(t,1); fn = ark(t,2);
        
        md1 = mean(data.obs.p1(st:fn,k)-data.obs.l1(st:fn,k),'omitnan');
        md2 = mean(data.obs.p2(st:fn,k)-data.obs.l2(st:fn,k),'omitnan');
        ml1 = mean(data.obs.l1(st:fn,k));
        ml2 = mean(data.obs.l2(st:fn,k));
        for i=st:fn
            data.obs.p1(i,k) = data.obs.l1(i,k) + md1 +...
                (2*f2^2*(1/(f1^2-f2^2)))*((data.obs.l1(i,k)-ml1)-(data.obs.l2(i,k)-ml2));
            data.obs.p2(i,k) = data.obs.l2(i,k) + md2 +...
                (2*f1^2*(1/(f1^2-f2^2)))*((data.obs.l1(i,k)-ml1)-(data.obs.l2(i,k)-ml2));
        end
    end
end
end