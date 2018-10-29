function [elip] = xyz2plh(cart,dopt)

if numel(cart)~= 3
    errdlg('Input matrice have to include maximum 3 component')
elseif size(cart,1) ~= 3
    cart = cart';
end

a  = 6378137.0;
f  = 1/298.257223563;
e2 = 2*f - f^2;

lam = atan2(cart(2),cart(1));
lam = mod(lam,2*pi);

p = sqrt(cart(1)^2 + cart(2)^2);

phi0 = atan(cart(3)/(p.*(1 - e2)));

while 1
    N   = a/sqrt(1 - (e2*(sin(phi0)^2)));
    h   = p/cos(phi0) - N;
    phi = atan((cart(3)/p)/(1 - (N/(N+h)*e2)));
    dphi = abs(phi - phi0);
    
    if dphi>10^-12
        phi0 = phi;
    else
        break
    end
end
if dopt==0
    elip = [phi;lam;h];
elseif dopt == 1
    t = (180/pi);
    elip = [(phi*t);(lam*t);h];
end

end

