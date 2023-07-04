%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NMSE performance vs Measurements
% SS-GAMP, OAMP-MMV, TS-OAMP-MA, SWOMP
% 2 bits quantization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
warning off;

%% Parameters Settings
para = parameter_init();

% lambda = para.c/para.fc;                                                  % wave length
% Dept = para.Nr/2*lambda;                                                  % antenna aperture
% Rayleigh_distance = 2*Dept^2/lambda;                                      % Rayleigh distance

%% Dictionnary Matrix
sin_index = -1:2/para.Nr:1-2/para.Nr;
D_far = exp(1i*pi*(0:para.Nr-1).' * sin_index)/sqrt(para.Nr);               % traditional dictionary

redunt = 2;
sin_index = -1:2/redunt/para.Nr_sub:1-2/redunt/para.Nr_sub;
D = exp(1i*pi*(0:para.Nr_sub-1).' * sin_index)/sqrt(para.Nr_sub*redunt);   
switch para.Num_sub
    case 8
        D_sub = blkdiag(D,D,D,D,D,D,D,D);
    case 6
        D_sub = blkdiag(D,D,D,D,D,D);  
    case 5
        D_sub = blkdiag(D,D,D,D,D); 
    case 4
        D_sub = blkdiag(D,D,D,D);
    case 3
        D_sub = blkdiag(D,D,D);     
    case 2
        D_sub = blkdiag(D,D);
    case 1
        D_sub = blkdiag(D);
end

%% Initialization
alg_type = 4;                                                               % number of algorithms
alg_set = [1:4];
var_set = [50:25:200];
epsilon_set = [0.06, 0.04, 0.025, 0.03, 0.035, 0.035, 0.037];               % for SWOMP
Nsim = [ones(1,length(var_set))*500];                                       % number of simulations
method = char('SS-GAMP', 'OAMP-MMV', 'TS-OAMP', 'SWOMP');

AER = zeros(alg_type,length(var_set));                                      % activity error rate
Pfa = zeros(alg_type,length(var_set));                                      % false alarm probability
Pmd = zeros(alg_type,length(var_set));                                      % miss detection probability
NMSE = zeros(alg_type,length(var_set),para.niter);                          % NMSE

for i_var=1:length(var_set)
    para.M = var_set(i_var);
    for i_sim = 1:Nsim(i_var)
        %% Channel Generation
        H = channel_generation_near_field_ULA_new(para);
        
        %% Activity Generation
        Index_active = randperm(para.K,para.Ka);                            % actice device set
        act_flag = zeros(para.K,1);
        act_flag(Index_active) = 1;                                         % activity indicators
        Mask = zeros(para.K,para.Nr);
        for i_k = 1:para.Ka
            Mask(Index_active(i_k),:) = 1;
        end
        
        %% Measurement Matrix
        F = dftmtx(para.K)/sqrt(para.K);
        S = F(randperm(para.K,para.M),:);                                   % pilot matrix
        
        %% Transmission
        H_space = H.*Mask;
        Z = S * H_space;
        noise = sqrt(para.Pn/2) * (randn(size(Z)) + 1i*randn(size(Z)));
        Y = Z + noise;                    
        para.adj =  sqrt(2)*0.675/3;                                        % to contorl the quantization range
        para.Ymax = para.adj*max(max(abs((Y)),[],'all'));
        Y_norm = Y / para.Ymax;                                             % normalization ( AGC)
        
        for alg_sel_ii = 1:length(alg_set)
            alg_sel = alg_set(alg_sel_ii);
            %% Quantization
            [Y_Q_R,~,~] = quantizer(real(Y_norm), para.Nbits, para.Nmax, -para.Nmax);                    % real part
            [Y_Q_I,delta_h, codebook_h] = quantizer(imag(Y_norm), para.Nbits, para.Nmax, -para.Nmax);    % iamginary part
            Y_Q = Y_Q_R+ 1i * Y_Q_I;                                        % complexy signal
            
            %% Add Mixed ADC
            para.Num_mix = para.Nr/ para.intvl;                             % total number of high-resoultion ADCs
            para.Num_mix_sub = para.Num_mix/para.Num_sub;
            mix_index = para.intvl/2:para.intvl:para.Nr;                    % uniformly distributed
            
            Y_wave = Y_Q;
            Y_wave(:,mix_index) = Y_norm(:,mix_index);                      % high-resolution ADCs with infinite quantization

            %% Joint Activity Detection and Channel Estimation
            para.Nvar2 = para.Pn/para.Ymax/para.Ymax;                       % normalized noise variance
            switch alg_sel
                case {1}  % SS-GAMP, sapce
                    [h_est, Xhat_final, rou, iterAMP, iterInAMP, v, iter_stop_h] = gamp_quantize_v2_search(Y_wave, S, para.Sa, codebook_h(end), codebook_h(1), para.Nbits, 10,...
                        para.niter/10, para.Nvar2, 0.5, 1e-5, 1e-4, 4, 1, 1, 1, H_space, para);
                    h_est = h_est * para.Ymax;
                    H_est = Xhat_final * para.Ymax;
                    h_est(:,:,end) = H_est;
                    
                    rou_mean = reshape(rou,para.K,[]);
                    act_hat = mean(rou_mean,2)>para.threshold;
                    
                case {2}  % OAMP-MMV, angle
                    [h_est, rou, ~, v_s_h5] = OAMP_MMV(Y_wave,S,codebook_h, para, D_far);
                    h_est = h_est * para.Ymax;
                    H_est = h_est(:,:,end)*D_far';
                    
                    oamp_th = 0.05*max(max(abs(H_est*D_far)));
                    act_hat = zeros(para.K,1);
                    act_hat(mean(abs(H_est),2) > oamp_th) = 1;
                    
                case {3}  % TS-OAMP, AD in space & CE in angle
                    [h_est, act_hat, ~, v_s_h5, para] = TS_OAMP(Y_wave,S,codebook_h, para, mix_index, D_sub);
                    h_est = h_est * para.Ymax;
                    H_est = h_est(:,:,para.niter);
                    
                case {4}   %  SWOMP, space
                    epsilon = epsilon_set(i_var);
                    [h_est, iter_num] = SW_OMP_Algorithm(Y_wave, S, epsilon);
                    h_est = h_est * para.Ymax ;
                    H_est = h_est;
                    h_est = repmat(h_est, 1, 1, para.niter);
                    
                    swomp_th = 0.01*max(max(abs(H_est)));
                    act_hat = zeros(para.K,1);
                    act_hat(mean(abs(H_est),2) > swomp_th) = 1;
                    
                otherwise
                    disp('error! inexistent algorithm');
            end
            
            %% AER
            AER_tmp = sum(abs(act_hat - act_flag),'all')/para.K;
            Pfa_tmp = sum(abs(act_hat(find(act_flag==0)) - act_flag(find(act_flag==0))),'all')/(para.K-para.Ka);
            Pmd_tmp = sum(abs(act_hat(find(act_flag==1)) - act_flag(find(act_flag==1))),'all')/para.Ka;
            AER(alg_sel,i_var) = AER(alg_sel,i_var) + AER_tmp;
            Pfa(alg_sel,i_var) = Pfa(alg_sel,i_var) + Pfa_tmp;
            Pmd(alg_sel,i_var) = Pmd(alg_sel,i_var) + Pmd_tmp;
            
            %% NMSE
            for i_nmse=1:para.niter
                switch alg_sel
                    case {1,4} % SS-GAMP, SWOMP
                        h_est_tmp = h_est(:,:,i_nmse);
                    case {2} % OAMP-MMV
                        h_est_tmp = h_est(:,:,i_nmse)*D_far';
                    case {3} % TS-OAMP-MA
                        h_est_tmp = h_est(:,:,i_nmse);
                end
                NMSE_temp = norm(h_est_tmp - reshape(H_space,para.K,[]),'fro')^2/norm(reshape(H_space,para.K,[]),'fro')^2;
                NMSE(alg_sel, i_var, i_nmse) = NMSE(alg_sel, i_var, i_nmse) + NMSE_temp;
            end
            NMSE_temp = norm(H_est - reshape(H_space,para.K,[]),'fro')^2/norm(reshape(H_space,para.K,[]),'fro')^2;
         
            %% display
            if mod(i_sim,min(10,Nsim(i_var))) == 0 %|| i_sim == 1
                fprintf('sim=%d, Pt=%d dBm, M=%d, Num of H-ADCs=%d, %s: AER=%.6f, AER_tmp=%.6f, NMSE=%.2f dB, NMSE_tmp=%.2f dB, Pfa=%.4f, Pmd=%.4f \n',...
                      i_sim, para.Ps_dBm, para.M, para.Nr/para.intvl, method(alg_sel,:),AER(alg_sel, i_var)/i_sim, AER_tmp,...
                      10*log10(NMSE(alg_sel, i_var, end)/i_sim),10*log10(NMSE_temp), Pfa(alg_sel, i_var)/i_sim, Pmd(alg_sel, i_var)/i_sim);
            end
        end
        
        if mod(i_sim, 200)==0
            Perform_tmp.NMSE = 10*log10(NMSE ./ repmat(i_sim,1,1,para.niter));
            Perform_tmp.AER = AER ./ i_sim;
            Perform_tmp.Pfa = Pfa ./ i_sim;
            Perform_tmp.Pmd = Pmd ./ i_sim;
            save Performance_tmp.mat Perform_tmp;
        end
    end
end
Perform.NMSE = 10*log10(NMSE ./ repmat(Nsim,1,1,para.niter));
Perform.AER = AER ./ Nsim;
Perform.Pfa = Pfa ./ Nsim;
Perform.Pmd = Pmd ./ Nsim;
disp('Finished all');

% save NMSEvsM_Bit2.mat Perform;
%% Plot
var_set = [50:25:200];
MarkerSize = 6;
LineWidth = 2;
Fontsize = 15;
figure
p2=plot(var_set,Perform.NMSE(1,:,end),'k-<','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p3=plot(var_set,Perform.NMSE(2,:,end),'g-*','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p4=plot(var_set,Perform.NMSE(3,:,end),'b-s','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p1=plot(var_set,Perform.NMSE(4,:,end),'m-d','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
set(gca,'FontSize',Fontsize);
gcf_set = [p1,p2,p3,p4];
legend(gcf_set,{'SWOMP','SS-GAMP','OAMP-MMV','TS-OAMP'},'interpreter','latex');
xlabel('Length of pilot sequence $M$','Fontsize',Fontsize,'interpreter','latex');
ylabel('NMSE','Fontsize',Fontsize,'interpreter','latex');
