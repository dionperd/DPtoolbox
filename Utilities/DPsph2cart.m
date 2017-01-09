function [x, y, z] = DPsph2cart(r,theta,phi,deg)

if deg
    theta = degtorad(theta);
    phi = degtorad(phi);
end
x = r.*cos(theta).*sin(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(phi);