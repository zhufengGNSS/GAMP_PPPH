function [data] = preprocess(data,options)

[data]  = elm_badclk(data);

[data] = decimation(data,options);

[data] = cal_sat(data);

[data] = elv_mask(data,options);

[data] = cs_detect(data,options);

if options.clkjump == 1
    [data] = clk_jmp2(data);
end

[data] = outlier(data);

if options.codsmth == 1
    [data] = smoothing(data);
end
end

