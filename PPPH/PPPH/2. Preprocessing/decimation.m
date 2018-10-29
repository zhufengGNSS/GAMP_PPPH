function [data] = decimation(data,options)

if (options.from<data.obs.ep(1,1)) || (options.to>data.obs.ep(end,1))
    errordlg('Observation file may not contain data between the chosen interval.','Epoch Interval Error');
    error   ('Observation file may not contain data between the chosen interval.');
end

f = (options.from - data.obs.ep(1,1));
l = (options.to - data.obs.ep(1,1));
fe = round(f/data.inf.time.int) + 1; %first epoch
le = round(l/data.inf.time.int) + 1; %last epoch

if le<size(data.obs.st,1)
    data.obs.p1((le+1):end,:) = [];
    data.obs.p2((le+1):end,:) = [];
    data.obs.l1((le+1):end,:) = [];
    data.obs.l2((le+1):end,:) = [];
    data.obs.ep((le+1):end,:) = [];
    data.obs.st((le+1):end,:) = [];
end

if fe~=1
    data.obs.p1(1:(fe-1),:) = [];
    data.obs.p2(1:(fe-1),:) = [];
    data.obs.l1(1:(fe-1),:) = [];
    data.obs.l2(1:(fe-1),:) = [];
    data.obs.ep(1:(fe-1),:) = [];
    data.obs.st(1:(fe-1),:) = [];
end
end

