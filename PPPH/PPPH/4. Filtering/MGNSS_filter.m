function [xs,kofs,pks] = MGNSS_filter(model,data,options)


[satlist] = dtr_satlist(data.obs);

[satno] = dtr_satno(data.obs);

rn = size(data.obs.st,1);

[freq,~] = frequencies;

[bp,~,sit] = dtr_sys(options);

xs  = zeros(bp+105,rn);
pks = zeros(bp+105,bp+105,rn);
kofs = zeros(5,5,rn);

n = 1;

for i=1:rn
    
    sls = satlist{i};
    
    sno = satno(i);
    
    pno = bp + sno;
    
    nk = n + (4*sno) - 1;
    meas = model(n:nk,:);
    
    n  = nk + 1;
    
    Q  = zeros(pno);
    
    if options.ProMod == 0
        Q(1,1) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
        Q(2,2) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
        Q(3,3) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
    end
    Q(4,4) = (options.NosClk*10^(options.NosClk2))*(data.inf.time.int);
    Q(5,5) = (options.NosTrop*10^(options.NosTrop2))*(data.inf.time.int);
    
    switch options.TroGrad
        case 0
            if bp==6
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==7
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(7,7) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==8
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(7,7) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            end
        case 1
            Q(6,6) = (1*10^(-12))*(data.inf.time.int);
            Q(7,7) = (1*10^(-12))*(data.inf.time.int);
            if bp==8
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==9
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(9,9) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==10
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(9,9) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(10,10) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            end
    end
    
    F = eye(pno);
    
    if i == 1 && options.InMethod == 0
        
        x1k = zeros(pno,1);
        
        p1k = zeros(pno);
        p1k(1,1) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(2,2) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(3,3) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(4,4) = (options.IntClk*10^(options.IntClk2))^2;
        p1k(5,5) = (options.IntTrop*10^(options.IntTrop2))^2;
        switch options.TroGrad
            case 0
                if bp==6
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==7
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(7,7) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==8
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(7,7) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                end
            case 1
                p1k(6,6) = (0*10^(1))^2;
                p1k(7,7) = (0*10^(1))^2;
                if bp==8
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==9
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(9,9) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==10
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(9,9) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(10,10) = (options.IntSTD*10^(options.IntSTD2))^2;
                end
        end
        
        for u=(bp+1):pno
            p1k(u,u) = (options.IntAmb*10^(options.IntAmb2))^2;
        end
    elseif i ~= 1
        
        x1k = zeros(pno,1);
        x1k(1:bp,1) = xs(1:bp,i-1);
        for k=1:sno
            snm = sls(k);
            x1k(bp+k,1) = xs(bp+snm,i-1);
        end
        x1k = F*x1k;
        
        p1k = zeros(pno);
        for r=1:size(p1k,1)
            for c=1:size(p1k,2)
                if r<(bp+1) && c<(bp+1)
                    p1k(r,c) = ps(r,c);
                elseif r<(bp+1) && c>bp
                    sn = sls(c-bp);
                    p1k(r,c) = ps(r,sn+bp);
                elseif r>bp && c<(bp+1)
                    sn = sls(r-bp);
                    p1k(r,c) = ps(sn+bp,c);
                else
                    f1 = sls(r-bp);
                    f2 = sls(c-bp);
                    p1k(r,c) = ps(f1+bp,f2+bp);
                end
            end
        end
        
        for k=1:sno
            snm = sls(k)+bp;
            if ps(snm,snm)==0
                p1k(k+bp,k+bp) = (options.IntAmb*10^(options.IntAmb2))^2;
            end
        end
        
        p1k = F*p1k*F' + Q;
    end
    
    if i == 1
        if options.InMethod == 1
            [ xk,pk,kof ] = least_sqr(sno,sls,pno,meas,sit,bp,freq,options);
        else
            % initial point
            if strcmp(options.ApMethod,'RINEX') % first choice is from RINEX
                x1k(1:3,1) = data.inf.rec.pos';
            elseif strcmp(options.ApMethod,'Specify')% second choice is from User specific
                x1k(1:3,1) = [options.AprioriX options.AprioriY options.AprioriZ];
            end
            [ xk,pk,kof ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
        end
    else
        [ xk,pk,kof ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
    end 
    
    kofs(:,:,i) = kof(:,:,1);
    
    xs(1:bp,i) = xk(1:bp,1);
    for k = 1:sno
       snm = sls(k) + bp;
       xs(snm,i) = xk(bp+k,1);
    end
    ps = zeros(bp + 105);
    for r=1:size(pk,1)
       for c=1:size(pk,2)
           if r<(bp+1) && c<(bp+1)
               ps(r,c) = pk(r,c);
           elseif r<(bp+1) && c>bp
               sn = sls(c-bp);
               ps(r,sn+bp) = pk(r,c);
           elseif r>bp && c<(bp+1)
               sn = sls(r-bp);
               ps(sn+bp,c) = pk(r,c);
           else
               sn1 = sls(r-bp); sn2 = sls(c-bp);
               ps(sn1+bp,sn2+bp) = pk(r,c);
           end
       end
    end
    pks(:,:,i) = ps;
end
xs(1:3,:) = xs(1:3,:);
end