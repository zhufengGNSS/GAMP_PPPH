function [rclk] = rel_clk(s_xyz,v_xyz)

c = 299792458; %m/s

rclk = ((-2)*((dot(s_xyz,v_xyz))/(c^2)))*c;

end

