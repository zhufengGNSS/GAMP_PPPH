function [xs] = ppph_complete(files, options)

[data] = data_hand(files,options);

[data] = preprocess(data,options);

[model] = nmodel(data,options);

[xs,~,~] = MGNSS_filter(model,data,options);

xs = xs';
end