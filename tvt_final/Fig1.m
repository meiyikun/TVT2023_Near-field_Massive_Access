clc;
clear;
close all;
rng(1);

%% Parameters Settings
Nr = 512;
K = 100;
fc = 6e9;
para.r_max = 100;
para.r_min = 10;
para.Ps_dBm = 5;
c = 3e8;
lambda = c/fc;
Dept = Nr/2*lambda;
r_max = para.r_max;
r_min = para.r_min;
p_ad = 2;
Ps_max = 10^(para.Ps_dBm/10)*1e-3;

sin_index = -1:2/Nr:1-2/Nr;
D_far = exp(1i*pi*(0:Nr-1).' * sin_index)/sqrt(Nr);

%% BS Location
location0 = (r_max-Dept)/2;
location_BS = [zeros(Nr,1), (location0:lambda/2:location0+lambda/2*(Nr-1)).'];

%% Devices Locations
location_UE_n = zeros(K,2);
theta = -pi*75/180+pi*150/180*rand(K,1);
radius = randi([r_min,r_max],K,1);
location_UE_n(:,1) = radius.*cos(theta);
location_UE_n(:,2) = r_max/2 + radius.*sin(theta);
% selected device
location_UE_n(25,:) = [48,52];

%% Rician Factor Setting
K_r = 10.^(1.3-0.003*radius);  
K_r_LoS = sqrt(K_r./(K_r+1));
K_r_NLoS = sqrt(1./(K_r+1));

%% Scatters
Lp_scatter = 3;                                                             % number of scatters
location_scatter_n = zeros(Lp_scatter, 2);
theta = -pi*75/180+pi*150/180*rand(Lp_scatter,1);
radius = randi([r_min,r_max],Lp_scatter,1);
location_scatter_n(:,1) = radius.*cos(theta);
location_scatter_n(:,2) = r_max/2 + radius.*sin(theta);
% selected scatters
location_scatter_n(1,:) = [25,16];
location_scatter_n(2,:) = [47,87];
location_scatter_n(3,:) = [40,27];

%% Channel
H_LoS_n = zeros(Nr, K);
H_NLoS_n = zeros(Nr, K);
pathloss_n = zeros(Nr, K);
Ps = zeros(K,1);
for i_UE = 1:K
    %% Large-scale fading
    BS_UE_vec_n = location_BS - location_UE_n(i_UE, :);
    d_BS_UE_n = sqrt(BS_UE_vec_n(:,1).^2+BS_UE_vec_n(:,2).^2);
    pathloss_n(:,i_UE) = (c./fc/4/pi./d_BS_UE_n).^2;
    
    %% Small-scale fading   
    % LoS  
    H_LoS_n(:, i_UE) = (exp(1i*2*pi*fc.*d_BS_UE_n/c));                   
    
    % NLos
    for i_Lp = 1:Lp_scatter
        BS_scatter_vec = location_BS - location_scatter_n(i_Lp, :);
        d_BS_scatter = sqrt(BS_scatter_vec (:,1).^2+BS_scatter_vec (:,2).^2);
        scatter_UE_vec = location_scatter_n(i_Lp, :) - location_UE_n(i_UE, :);
        d_scatter_UE = sqrt(scatter_UE_vec(:,1).^2+scatter_UE_vec(:,2).^2);
        H_NLoS_n(:, i_UE) = H_NLoS_n(:, i_UE) + exp(1i*2*pi*fc*(d_scatter_UE+d_BS_scatter)/c);
    end
    
    %% Transmission Power    
    Ps(i_UE) = Ps_max/r_max^p_ad*d_BS_UE_n(Nr/2).^p_ad;
end
H_LoS_n = K_r_LoS.*H_LoS_n.';
H_NLoS_n = K_r_NLoS.*H_NLoS_n.';
H_n = sqrt(Ps.*pathloss_n.').*(H_LoS_n + H_NLoS_n);

Index_active = 25;
Mask = zeros(K,Nr);
Mask(Index_active,:) = 1;
H_plot_n = H_n.*Mask*D_far;
H_plot_n = H_plot_n/max(abs(H_plot_n),[],'all');

%% Far-Field Channel
beta=21;
location_scatter_f = zeros(Lp_scatter, 2);
location_UE_f = [48*beta,52*beta-50*beta+50];
location_scatter_f(1,:) = [25*beta,16*beta-50*beta+50];
location_scatter_f(2,:) = [47*beta,87*beta-50*beta+50];
location_scatter_f(3,:) = [40*beta,27*beta-50*beta+50];

BS_UE_vec_f = location_BS - location_UE_f;
d_BS_UE_f = sqrt(BS_UE_vec_f(:,1).^2+BS_UE_vec_f(:,2).^2);
pathloss_f = (c./fc/4/pi./d_BS_UE_f).^2;
H_LoS_f = exp(1i*2*pi*fc.*d_BS_UE_f/c);                   

H_NLoS_f = zeros(Nr,1);
for i_Lp = 1:Lp_scatter
    BS_scatter_vec_f = location_BS - location_scatter_f(i_Lp, :);
    d_BS_scatter_f = sqrt(BS_scatter_vec_f (:,1).^2+BS_scatter_vec_f (:,2).^2);
    scatter_UE_vec_f = location_scatter_f(i_Lp, :) - location_UE_f;
    d_scatter_UE_f = sqrt(scatter_UE_vec_f(:,1).^2+scatter_UE_vec_f(:,2).^2);
    H_NLoS_f= H_NLoS_f + exp(1i*2*pi*fc*(d_scatter_UE_f+d_BS_scatter_f)/c);
end 
Ps = Ps_max/r_max^p_ad*d_BS_UE_f(Nr/2).^p_ad;

H_LoS_f = K_r_LoS(Index_active).*H_LoS_f.';
H_NLoS_f = K_r_NLoS(Index_active).*H_NLoS_f.';
H_f = sqrt(Ps.*pathloss_f.').*(H_LoS_f + H_NLoS_f);

H_plot_f = H_f*D_far;
H_plot_f = H_plot_f/max(abs(H_plot_f),[],'all');

%% plot
figure;
plot(sin_index,abs(H_plot_n(Index_active,:))); hold on;
plot(sin_index,abs(H_plot_f));
Fontsize = 15;
set(gca,'xtick',-1:0.2:1,'FontSize',Fontsize);
xlabel('sin($\theta$)','Fontsize',Fontsize,'interpreter','latex');
ylabel('Normalized Amplitude','Fontsize',Fontsize,'interpreter','latex');
legend('Near-field','Far-field','interpreter','latex');



