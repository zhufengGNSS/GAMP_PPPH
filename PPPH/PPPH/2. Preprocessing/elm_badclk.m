function [data] = elm_badclk(data)

intr = data.inf.time.clkint;
for i=1:size(data.clk,2)
    if (any(~isnan(data.clk(:,i)))) && (any(isnan(data.clk(:,i)) | data.clk(:,i)==999999.999999))
        loc = find(isnan(data.clk(:,i)));
        for t=loc'
            sod = (t-1)*intr;
            ep  = data.obs.ep(:,1)==sod;
            data.obs.st(ep,i) = 0;
        end
    elseif any(isnan(data.clk(:,i)))
        data.obs.st(:,i) = 0;
    end
end

end

