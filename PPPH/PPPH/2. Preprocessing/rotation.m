function [out] = rotation(position,angle,axis)

narginchk(3,3)

if numel(position) ~= 3
     error('Matrix dimension shold be 3xN or Nx3.')
end


if numel(angle) ~= 1 || numel(axis) ~= 1
    error('Angle and axis should be scalar.')
end


switch axis
    case 1
        rot = [1 0 0; 0 cosd(angle) sind(angle); 0 -sind(angle) cosd(angle)];
    case 2
        rot = [cosd(angle) 0 -sind(angle); 0 1 0; sind(angle) 0 cosd(angle)];
    case 3 
        rot = [cosd(angle) sind(angle) 0; -sind(angle) cosd(angle) 0; 0 0 1];
end


[r,~] = size(position);

switch r
    case 1
        xout = rot*(position');
    case 3
        xout = rot*(position);
end

out = xout';
nargoutchk(1,1)
end

