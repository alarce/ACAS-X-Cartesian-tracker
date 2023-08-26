function [trk] = InitializeCartesianTracker(trk, r_slant, Chi_abs, z_rel)
%trk: modeSTrackFile pe8
global params;
R_k = params.surveillance.horizontal.cartfilter.R; %noise matrix for measurement
U = params.surveillance.horizontal.cartfilter.U; %
kappa = params.surveillance.horizontal.cartfilter.kappa; %initial covariance inflation, applied to obtain positive semidefinite matrixes
gamma = params.surveillance.horizontal.cartfilter.gamma; %covariance inflation factor

%logic for no bearing track file -> NA

if (abs(z_rel) > abs(r_slant))
    %para asegurar que todas las matrices son positivas semi definidas
    r_slant = abs(z_rel); 
end
[O_pol, w] = SigmaPointSample([r_slant; Chi_abs], R_k, 3/2);%observed en polares
O_xy = zeros(2, 5); %o:observed en cartesianas x e y
zero_ground_range = false;
for i = 1:5 %5 sigma points
    if (abs(O_pol(1,i)) <= abs(z_rel))
    zero_ground_range = true;   %comprobar que para cada uno de los puntos sigma
    end
    O_xy(:,i) = ConvertToCartesian(O_pol(1,i), O_pol(2,i), z_rel);
end
[mu_xy, Sigma_xy] = WeightedMeanAndCovariance(O_xy,w);
 %mu_xy = ConvertToCartesian(r_slant, Chi_abs, z_rel)  %importante añadido por Alvaro
if (zero_ground_range)
    Sigma_xy = Sigma_xy + gamma*eye(2);
end

%this is not acas code, initialize with true values
global trueXInitial;
global trueYInitial;
global trueXVelInitial;
global trueYVelInitial;
mu_xy = [trueXInitial; trueYInitial];
%end of not acas code

%trk.mu_cart = [mu_xy; trueXVelInitial; trueYVelInitial]; %ALVARO not acas code
trk.mu_cart = [mu_xy; zeros(2,1)]; %x,y,vx,vy %true acas code
%trk.mu_cart = [mu_xy; -800; 0] %alvaro
%Sigma_xy(1,2) = 0%alvaro
%Sigma_xy(2,1) = 0%alvaro
trk.Sigma_cart = block_diag(Sigma_xy, U) + kappa * eye(4);
trk.updates_cart = 1;
end



    
    
    