function [mu, Sigma] = UpdateCartesianTracker(r_slant, Chi_abs, z_rel, mu_s, Sigma_s)

global params;
R_k = params.surveillance.horizontal.cartfilter.R; %measurement noise
%[mu_aug, Sigma_aug] = [[mu_s; zeros(2)], block_diag(Sigma_s,R_k)];
mu_aug = [mu_s; zeros(2,1)]; %augmented
Sigma_aug = block_diag(Sigma_s,R_k);

[S_aug, w_aug] = SigmaPointSample(mu_aug, Sigma_aug, 9/2);
gamma = zeros(2,13);
for i = 1:13 %13 (2m+1=2*6+1) sigma points per state variable
    gamma(:,i) = AugmentedStateToObservation(S_aug(:,i), z_rel);
end
[mu_o,Sigma_o] = WeightedMeanAndCovariance(gamma,w_aug);
K = zeros(length(mu_s), length(mu_o));
for i = 1:length(w_aug)
    K = K + w_aug(i)*(S_aug(1:4,i) - mu_s)*(gamma(:,i) - mu_o)';
end
K = K * inv(Sigma_o);
o = ConvertToCartesian(r_slant, Chi_abs, z_rel);
mu = mu_s + K*(o - mu_o);
Sigma = Sigma_s - K*Sigma_o*K';
Sigma = (Sigma + Sigma') ./ 2;

