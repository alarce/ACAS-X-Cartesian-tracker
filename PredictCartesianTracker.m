function [mu_s, Sigma_s] = PredictCartesianTracker(mu, Sigma, dt)

global params;
min_extrap_toa_step = params.surveillance.min_extrap_toa_step;
Q = params.surveillance.horizontal.cartfilter.Q; %process noise %512
if (dt < min_extrap_toa_step)
    mu_s = mu;
    Sigma_s = Sigma;
else
    F = [1 0 dt 0; ...
    0 1 0 dt; ...
    0 0 1 0; ...
    0 0 0 1];
    Gdt = 0.5*dt^2;
    G = [Gdt 0 ; ...
    0 Gdt; ...
    dt 0 ; ...
    0 dt ];
    [mu_s, Sigma_s] = PredictKalmanFilter(mu, Sigma, F, G, Q);
    
    %debug
    acasx_x_predicted = mu_s(1); %debug
    acasx_y_predicted = mu_s(2); %debug
end

