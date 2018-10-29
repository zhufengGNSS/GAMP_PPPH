function [bp,ap,sit] = dtr_sys(options)

if options.TroGrad == 0
    base = 5;
elseif options.TroGrad == 1
    base = 7;
end

if options.system.gps == 1
    if options.system.glo == 1
        if options.system.gal == 1
            if options.system.bds == 1 
                bp  =  base + 3;
                ap  = 92;
                sit = 12;
            else                      
                bp  =  base + 2;
                ap  = 78;
                sit =  8;
            end
        elseif options.system.bds == 1 
            bp  =  base + 2;
            ap  = 72;
            sit =  9;
        else                          
            bp  =  base + 1;
            ap  = 58;
            sit =  3;
        end
    elseif options.system.gal == 1
        if options.system.bds == 1     
            bp  =  base + 2;
            ap  = 66;
            sit = 10;
        else                           
            bp  =  base + 1;
            ap  = 52;
            sit =  4;
        end
    elseif options.system.bds == 1     
        bp  =  base + 1;
        ap  = 46;
        sit =  5;
    else                               
        bp  =  base;
        ap  = 32;
        sit =  1;
    end
elseif options.system.glo == 1
    if options.system.gal == 1
        if options.system.bds == 1     
            bp  =  base + 2;
            ap  = 60;
            sit = 11;
        else                           
            bp  =  base + 1;
            ap  = 46;
            sit =  6;
        end
    elseif options.system.bds == 1    
        bp  =  base + 1;
        ap  = 40;
        sit =  7;
    else                               
        bp  =  base;
        ap  = 26;
        sit =  2;
    end
else
    errordlg('Process must include GPS or GLONASS satellites','Navigation System Error')
end

end

