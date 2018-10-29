function [arc] = arc_dtr(obs)

sn  = size(obs.st,2);
arc = cell(1,sn);

for k=1:sn
    
    dc = find(isnan(obs.p1(:,k)) & obs.st(:,k)==1);
    if any(dc)
        obs.st(dc,k) = 0;
    end
    
    row = find(obs.st(:,k)==1);
    
    brk = find(diff(row)~= 1);
    
    frs = [1;(brk)+1]; lst = [brk;size(row,1)];
    als = [frs lst];
    
    if size((als),1)>1
        for t=size((als),1):-1:1
            if (als(t,2)-als(t,1))<10
                als(t,:) = [];
            end
        end
    else
        if (als(1,2)-als(1,1))<10
            als(1,:) = [];
        end
    end
    
    arn = NaN(size(als,1),2);
    for i=1:size(als,1)
        arn(i,1) = row(als(i,1)); arn(i,2) = row(als(i,2));
    end
    arc{k} = arn;
end

end

