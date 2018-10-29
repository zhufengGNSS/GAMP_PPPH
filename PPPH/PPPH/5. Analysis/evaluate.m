function [n,e,u,CT,thrD,rms] = evaluate(xs,ref)
    set = xs(1:3,:) - ref;
    [elip] = xyz2plh(ref,1);
    lat = elip(1);
    lon = elip(2);

    % [lat,lon,~]=xyz2plh(ref(1),ref(2),ref(3));
    % Transformation matrix from Global to Local
    A = [-sind(lat)*cosd(lon), -sind(lat)*sind(lon), cosd(lat);...
         -sind(lon),            cosd(lon),           0;...
          cosd(lat)*cosd(lon),  cosd(lat)*sind(lon), sind(lat)];
        
%     A = [(-sind(lat)*cosd(lon)) (-sind(lon)) (cosd(lat)*cosd(lon));...
%          (-sind(lat)*sind(lon)) ( cosd(lon)) (cosd(lat)*sind(lon));...
%          (           cosd(lat)) (         0) (          sind(lat))];

    loc = zeros(size(set,2),3);
    for i = 1:size(set,2)
        loc(i,:) = (A*set(:,i))';
    end

    n = loc(:,1); e = loc(:,2); u = loc(:,3);
    thrD = sqrt(n.^2 + e.^2 + u.^2);
    CT = 1;
    dur = 19;
    for i=1:size(set,2)
        if i+dur<size(thrD,1)
            val = thrD(i:i+dur,1)>0.1;
            if sum(val)==0
                CT = i;
                break
            end
        else
            val = thrD(i:end,1)>0.1;
            if sum(val)==0
                CT = i;
                break
            end
        end
    end
    rms = zeros(3,1);
    rms(1,1) = sqrt(mean(n(CT:end,1).^2));
    rms(2,1) = sqrt(mean(e(CT:end,1).^2));
    rms(3,1) = sqrt(mean(u(CT:end,1).^2));
end