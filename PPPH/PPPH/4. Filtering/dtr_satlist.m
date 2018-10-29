function [satlist] = dtr_satlist(obs)

m = size(obs.st,1);
satlist = cell(m,1);

for i=1:m
    satlist{i,1} = find(obs.st(i,:)==1);
end

end
