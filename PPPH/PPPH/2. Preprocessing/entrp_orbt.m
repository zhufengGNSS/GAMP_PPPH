function [out1,out2] = entrp_orbt(nep,gap,dat)

n = (nep/gap) + 6;

st = floor(n) - 4;
fn = ceil(n) + 4;
kern = dat(st:fn,1);
nt   = rem(n,1) + 5;
out1 = lag(nt,kern);
out2 = velo(nt,kern,gap);

end