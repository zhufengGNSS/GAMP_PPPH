function [satno] = dtr_satno(obs)

m = size(obs.st,1);
satno = zeros(m,1);

for i=1:m
    satno(i,1) = sum(obs.st(i,:));
end

end

