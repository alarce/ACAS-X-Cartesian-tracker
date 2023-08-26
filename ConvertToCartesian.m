function [cart] = ConvertToCartesian(r_slant, Chi_abs, z_rel)

rho = GroundRange(r_slant, z_rel);

cart = [rho*sin(Chi_abs); rho*cos(Chi_abs)];