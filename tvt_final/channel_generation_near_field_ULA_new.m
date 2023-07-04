function [H] = channel_generation_near_field_ULA_new(para)
%% Parameters Settings
Nr = para.Nr;
K = para.K;
fc = para.fc;
c = 3e8;
lambda = c/fc;
Dept = Nr/2*lambda;
r_max = para.r_max;
r_min = para.r_min;
p_ad = 2;
Ps_max = 10^(para.Ps_dBm/10)*1e-3;

%% BS Location
location0 = (r_max-Dept)/2;
location_BS = [zeros(Nr,1), (location0:lambda/2:location0+lambda/2*(Nr-1)).'];

%% Devices Locations
location_UE = zeros(K,2);
theta = -pi*75/180 + pi*150/180*rand(K,1);
radius = randi([r_min,r_max],K,1);
location_UE(:,1) = radius.*cos(theta);
location_UE(:,2) = r_max/2 + radius.*sin(theta);

%% Rician Factor Setting
K_r = 10.^(1.3-0.003*radius);  
K_r_LoS = sqrt(K_r./(K_r+1));
K_r_NLoS = sqrt(1./(K_r+1));

%% Scatters
Lp_scatter = 5;                                                             % number of scatters
location_scatter = zeros(Lp_scatter, 2);
theta = -pi*75/180+pi*150/180*rand(Lp_scatter,1);
radius = randi([r_min,r_max],Lp_scatter,1);
location_scatter(:,1) = radius.*cos(theta);
location_scatter(:,2) = r_max/2 + radius.*sin(theta);

%% Channel
H_LoS = zeros(Nr, K);
H_NLoS = zeros(Nr, K);
pathloss = zeros(Nr, K);
Ps = zeros(K,1);
for i_UE = 1:K
    %% Large-scale fading
    BS_UE_vec = location_BS - location_UE(i_UE, :);
    d_BS_UE = sqrt(BS_UE_vec(:,1).^2+BS_UE_vec(:,2).^2);
    pathloss(:,i_UE) = (c./fc/4/pi./d_BS_UE).^2;
    
    %% Small-scale fading   
    % LoS  
    H_LoS(:, i_UE) = (exp(1i*2*pi*fc.*d_BS_UE/c));                   
    
    % NLos
    for i_Lp = 1:Lp_scatter
        BS_scatter_vec = location_BS - location_scatter(i_Lp, :);
        d_BS_scatter = sqrt(BS_scatter_vec (:,1).^2+BS_scatter_vec (:,2).^2);
        scatter_UE_vec = location_scatter(i_Lp, :) - location_UE(i_UE, :);
        d_scatter_UE = sqrt(scatter_UE_vec(:,1).^2+scatter_UE_vec(:,2).^2);
        H_NLoS(:, i_UE) = H_NLoS(:, i_UE) + exp(1i*2*pi*fc*(d_scatter_UE+d_BS_scatter)/c);
    end
    
    %% Transmission Power    
    Ps(i_UE) = Ps_max/r_max^p_ad*d_BS_UE(Nr/2).^p_ad;
end
H_LoS = K_r_LoS.*H_LoS.';
H_NLoS = K_r_NLoS.*H_NLoS.';
H = sqrt(Ps.*pathloss.').*(H_LoS + H_NLoS);
