function [out1,out2] = entrp(nep,gap,dat)

min = 1;
max = 86400/gap;
n = (nep/gap) + 1;

if n<=(min+5)
    kern = dat(min:min+9,1);
    nt   = rem(n,(min+5));
    out1 = lag( nt,kern);
    out2 = velo(nt,kern,gap);
elseif n>(max-5)
    kern = dat(max-9:max,1);
    nt   = rem(n,(max-5)) + 5;
    out1 = lag(nt,kern);
    out2 = velo(nt,kern,gap);
else
    st = floor(n) - 4;
    fn = ceil(n) + 4;
    kern = dat(st:fn,1);
    nt   = rem(n,1) + 5;
    out1 = lag(nt,kern);
    out2 = velo(nt,kern,gap);
end
end