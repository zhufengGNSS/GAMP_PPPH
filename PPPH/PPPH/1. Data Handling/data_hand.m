function [data] = data_hand(files,options)

narginchk(2,2)

[obs,inf] = read_obsf(files,options);

[sat,inf] = r_sp3(files.orbit,options,inf);

options.entrp = 0; 
if ~isempty(files.orbitb) && ~isempty(files.orbitb)
    
    [satb,~] = r_sp3(files.orbitb,options,inf);
    
    [sata,~] = r_sp3(files.orbita,options,inf);

    if (size(satb.sp3,1)==size(sata.sp3,1))&&(size(satb.sp3,1)==size(sat.sp3,1))...
            &&(size(sat.sp3,1)==size(sata.sp3,1))
        sat.sp3 = vertcat(satb.sp3(end-4:end,:,:),...
                          sat.sp3(:,:,:),...
                          sata.sp3(1:5,:,:));
    end
    options.entrp = 1;
end

if strcmp(options.clock,'Clock File')
    [clk] = r_clck(files.clock,options);
    inf.time.clkint = options.clck_int;
elseif strcmp(options.clock,'Sp3 File')
    [clk] = sat.sp3(:,4,:);
    inf.time.clkint = inf.time.sp3int;
end

[atx] = r_antx(files.anten,inf,options);

data.inf  = inf;
data.robs = obs;
data.obs  = obs;
data.sat  = sat;
data.clk  = clk;
data.atx  = atx;
data.opt  = options;
data.files = files;
end