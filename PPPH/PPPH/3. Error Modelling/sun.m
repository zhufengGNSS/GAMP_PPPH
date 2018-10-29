function [s_pos] = sun(mjd)

AU  = 149597870700;                                                   % met
d2r = pi/180;
fday = mjd - floor(mjd);
JDN  = mjd - 15019.5;

v1   = mod((279.696678 + 0.9856473354*JDN),360);                      % deg
gstr = mod((279.690983 + 0.9856473354*JDN + 360*fday + 180),360);     % deg
g    = (mod((358.475845 + 0.9856002670*JDN),360))*d2r;                % rad

slong = v1 + (1.91946 - 0.004789*JDN/36525)*sin(g)...
      + 0.020094*sin(2*g);                                            % deg
obliq = (23.45229 - 0.0130125*JDN/36525)*d2r;                         % rad

slp  = (slong - 0.005686)*d2r;                                        % rad
snd  = sin(obliq)*sin(slp);
csd  = sqrt(1 - snd^2);
sdec = atan2d(snd,csd);                                               % deg

sra = 180 - atan2d((snd/csd/tan(obliq)),(-cos(slp)/csd));             % deg

s_pos = [(cosd(sdec)*cosd(sra)*AU);...
         (cosd(sdec)*sind(sra)*AU);...
         (sind(sdec)*AU)];

s_pos = rotation(s_pos,gstr,3);                                       % met

end

