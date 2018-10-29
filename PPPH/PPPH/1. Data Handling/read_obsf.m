function [ obs,inf ] = read_obsf(files,options)

[ver] = r_rnxvers(files.rinex);

if ver>=3
    
    [inf] = r_rnxheadv3(files.rinex);
    
    if options.dcb == 1
        [dcb] = r_dcb(files.dcb);
        
        [obs] = r_rnxobsv3(files.rinex,inf,options,dcb);
    else
        [obs] = r_rnxobsv3(files.rinex,inf,options);
    end
elseif ver>=2
    
    [inf] = r_rnxheadv2(files.rinex);
    
    if options.dcb == 1
        [dcb] = r_dcb(files.dcb);
        
        [obs] = r_rnxobsv2(files.rinex,inf,options,dcb);
    else
        
        [obs] = r_rnxobsv2(files.rinex,inf,options);
    end
else
    errordlg('RINEX version is not valid !','RINEX version error');
    error('RINEX version is not valid !');
end
end

