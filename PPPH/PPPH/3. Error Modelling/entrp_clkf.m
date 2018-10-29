function [out1] = entrp_clkf(nep,gap,dat)

min = 1;
max = 86400/gap;
n = (nep/gap) + 1;
ns = (10*60)/gap;
nt = (20*60)/gap;

if (n-ns)<=min
    st = min; fn=min+(nt-1);
    t  = st:fn;
    kern = dat(st:fn,1);
    [p,~,mu] = polyfit(t', kern, 2);
    out1 = polyval(p,n,[],mu);
elseif (n+ns)>max
    st = max-(nt-1); fn=max;
    t  = st:fn;
    kern = dat(st:fn,1);
    [p,~,mu] = polyfit(t', kern, 2);
    out1 = polyval(p,n,[],mu);
else
    st = floor(n) - (ns-1);
    fn = ceil(n)  + (ns-1);
    t  = st:fn;
    kern = dat(st:fn,1);
    [p,~,mu] = polyfit(t', kern, 2);
    out1 = polyval(p,n,[],mu);
end
end