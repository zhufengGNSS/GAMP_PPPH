function [p,T,dT,e,ah,aw,undu] = gpt2 (dmjd,dlat,dlon,hell,nstat,it)

dmjd1 = dmjd-51544.5;

gm = 9.80665;

dMtr = 28.965*10^-3;

Rg = 8.3143;

if (it==1)
    cosfy = 0;
    coshy = 0;
    sinfy = 0;
    sinhy = 0;
else 
    cosfy = cos(dmjd1/365.25*2*pi);
    coshy = cos(dmjd1/365.25*4*pi);
    sinfy = sin(dmjd1/365.25*2*pi);
    sinhy = sin(dmjd1/365.25*4*pi);
end

fid = fopen('gpt2_5.grd','r');


line = fgetl(fid);

for n = 1:2592
    
    line = fgetl(fid);
    vec = sscanf(line,'%f');
    
    pgrid(n,1:5)  = vec(3:7);          
    Tgrid(n,1:5)  = vec(8:12);         
    Qgrid(n,1:5)  = vec(13:17)./1000;  
    dTgrid(n,1:5) = vec(18:22)./1000;  
    u(n)          = vec(23);           
    Hs(n)         = vec(24);           
    ahgrid(n,1:5) = vec(25:29)./1000;  
    awgrid(n,1:5) = vec(30:34)./1000;  
    
end
fclose (fid);

for k = 1:nstat
    
    if dlon(k) < 0
        plon = (dlon(k) + 2*pi)*180/pi;
    else
        plon = dlon(k)*180/pi;
    end
    
    ppod = (-dlat(k) + pi/2)*180/pi;
    
    ipod = floor((ppod+5)/5); 
    ilon = floor((plon+5)/5);
    
    diffpod = (ppod - (ipod*5 - 2.5))/5;
    difflon = (plon - (ilon*5 - 2.5))/5;
    
    if ipod == 37
        ipod = 36;
    end
    
    indx(1) = (ipod - 1)*72 + ilon;
    
    bilinear = 0;
    if ppod > 2.5 && ppod < 177.5 
           bilinear = 1;          
    end          
    
    if bilinear == 0

        ix = indx(1);
        
        undu(k) = u(ix);
        hgt = hell(k)-undu(k);
        
        T0 = Tgrid(ix,1) + ...
             Tgrid(ix,2)*cosfy + Tgrid(ix,3)*sinfy + ...
             Tgrid(ix,4)*coshy + Tgrid(ix,5)*sinhy;
        p0 = pgrid(ix,1) + ...
             pgrid(ix,2)*cosfy + pgrid(ix,3)*sinfy+ ...
             pgrid(ix,4)*coshy + pgrid(ix,5)*sinhy;
         
        Q = Qgrid(ix,1) + ...
            Qgrid(ix,2)*cosfy + Qgrid(ix,3)*sinfy+ ...
            Qgrid(ix,4)*coshy + Qgrid(ix,5)*sinhy;
        
        dT(k) = dTgrid(ix,1) + ...
                dTgrid(ix,2)*cosfy + dTgrid(ix,3)*sinfy+ ...
                dTgrid(ix,4)*coshy + dTgrid(ix,5)*sinhy;
            
        redh = hgt - Hs(ix);

        
        T(k) = T0 + dT(k)*redh - 273.15;
        
        dT(k) = dT(k)*1000;
        
        Tv = T0*(1+0.6077*Q);
        
        c = gm*dMtr/(Rg*Tv);
        
        p(k) = (p0*exp(-c*redh))/100;
        
        e(k) = (Q*p(k))/(0.622+0.378*Q);
        
        ah(k) = ahgrid(ix,1) + ...
                ahgrid(ix,2)*cosfy + ahgrid(ix,3)*sinfy+ ...
                ahgrid(ix,4)*coshy + ahgrid(ix,5)*sinhy;
            
        aw(k) = awgrid(ix,1) + ...
                awgrid(ix,2)*cosfy + awgrid(ix,3)*sinfy + ...
                awgrid(ix,4)*coshy + awgrid(ix,5)*sinhy;           
                    
     else
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
        if ilon1 == 73
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 72;
        end
        
        indx(2) = (ipod1 - 1)*72 + ilon;  
        indx(3) = (ipod  - 1)*72 + ilon1; 
        indx(4) = (ipod1 - 1)*72 + ilon1;
        
        for l = 1:4
            
            undul(l) = u(indx(l));
            hgt = hell(k)-undul(l);
            
            T0 = Tgrid(indx(l),1) + ...
                 Tgrid(indx(l),2)*cosfy + Tgrid(indx(l),3)*sinfy + ...
                 Tgrid(indx(l),4)*coshy + Tgrid(indx(l),5)*sinhy;
            p0 = pgrid(indx(l),1) + ...
                 pgrid(indx(l),2)*cosfy + pgrid(indx(l),3)*sinfy + ...
                 pgrid(indx(l),4)*coshy + pgrid(indx(l),5)*sinhy;
             
            Ql(l) = Qgrid(indx(l),1) + ...
                    Qgrid(indx(l),2)*cosfy + Qgrid(indx(l),3)*sinfy + ...
                    Qgrid(indx(l),4)*coshy + Qgrid(indx(l),5)*sinhy;
                
            Hs1 = Hs(indx(l));
            redh = hgt - Hs1;
            
            dTl(l) = dTgrid(indx(l),1) + ...
                     dTgrid(indx(l),2)*cosfy + dTgrid(indx(l),3)*sinfy + ...
                     dTgrid(indx(l),4)*coshy + dTgrid(indx(l),5)*sinhy;
                 
            Tl(l) = T0 + dTl(l)*redh - 273.15;
            
            Tv = T0*(1+0.6077*Ql(l));  
            c = gm*dMtr/(Rg*Tv);
            
            pl(l) = (p0*exp(-c*redh))/100;
            
            ahl(l) = ahgrid(indx(l),1) + ...
                     ahgrid(indx(l),2)*cosfy + ahgrid(indx(l),3)*sinfy + ...
                     ahgrid(indx(l),4)*coshy + ahgrid(indx(l),5)*sinhy;
                 
            awl(l) = awgrid(indx(l),1) + ...
                     awgrid(indx(l),2)*cosfy + awgrid(indx(l),3)*sinfy + ...
                     awgrid(indx(l),4)*coshy + awgrid(indx(l),5)*sinhy;
            
        end
            
        dnpod1 = abs(diffpod); 
        dnpod2 = 1 - dnpod1;   
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
        
        R1 = dnpod2*pl(1)+dnpod1*pl(2);
        R2 = dnpod2*pl(3)+dnpod1*pl(4);
        p(k) = dnlon2*R1+dnlon1*R2;
        
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2);
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4);
        T(k) = dnlon2*R1+dnlon1*R2;
        
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
        
        R1 = dnpod2*Ql(1)+dnpod1*Ql(2);
        R2 = dnpod2*Ql(3)+dnpod1*Ql(4);
        Q = dnlon2*R1+dnlon1*R2;
        e(k) = (Q*p(k))/(0.622+0.378*Q);
        
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
        
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
                    
    end 
end

  
