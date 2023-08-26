function [trk] = AdvanceCartesianTrack(trk, r_slant, Chi_abs, z_rel, toa)

global params;
alt_threshold = params.surveillance.horizontal.cartfilter.altitude_threshold;
max_outliers_normal = params.surveillance.horizontal.cartfilter.max_outlier_detections;
max_outliers_reduced = params.surveillance.horizontal.cartfilter.max_outlier_detections_reduced;
outlier_thresh = 90; %params.surveillance.horizontal.cartfilter.outlier_threshold;
detections_to_recover = params.surveillance.horizontal.cartfilter.detections_to_recover;
R_q = params.surveillance.horizontal.cartfilter.R;
I_xy = [1,2];
max_outliers = max_outliers_normal;

dt = toa - trk.toa_cart;
[mu_s, Sigma_s] = PredictCartesianTracker(trk.mu_cart, trk.Sigma_cart, dt); %debug uncomment
[mu_pol, Sigma_pol] = LinearTransform(mu_s(I_xy), Sigma_s(I_xy,I_xy));
Chi_abs = WrapToPi(Chi_abs);
o_pol = [GroundRange(r_slant,z_rel); Chi_abs];
d_pol = [o_pol(1) - mu_pol(1); AngleDifference(mu_pol(2), o_pol(2))]; %residuo

%if (isnan(r_slant)) || (isnan(Chi_abs)) || (abs(z_rel) - r_slant >= alt_threshold) || (IsOutlier(d_pol, 0.0, Sigma_pol + R_q, outlier_thresh)) || (trk.updates_cart == 0)
%    if (trk.odc_cart > 0) && (trk.updates_cart > 1) && (trk.valid_cart)
%        trk.odc_cart = trk.odc_cart - 1;
%    else
%        if (trk.updates_cart==1)
%            trk.valid_cart = false;
%        end
%        InitializeCartesianTracker(trk,r_slant,Chi_abs,z_rel);
%        trk.odc_cart = max_outliers;
%        trk.toa_cart = toa;
%    end
%else
    [trk.mu_cart, trk.Sigma_cart] = UpdateCartesianTracker(r_slant, Chi_abs, z_rel, mu_s, Sigma_s); %debug uncomment
    %trk.mu_cart  %uncomment to debug
    %trk.Sigma_cart %uncomment to debug
    trk.odc_cart = max_outliers;
    trk.updates_cart = trk.updates_cart + 1;
    if (trk.updates_cart >= detections_to_recover)
        trk.valid_cart = true;
    end
    trk.toa_cart = toa;
%end

    
    
    
    