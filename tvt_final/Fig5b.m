%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NMSE performance vs Power Pt
% TS-OAMP w/o sub-array, TS-OAMP-MA w/ sub-array
% 2 bits quantization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
warning off;

%% Parameters Settings
para = parameter_init();
  
% lambda = para.c/para.fc;                                                    % wave length
% Dept = para.Nr/2*lambda;                                                    % antenna aperture
% Rayleigh_distance = 2*Dept^2/lambda;                                        % Rayleigh distance

%% Initialization
alg_type = 2;                                                               % number of algorithms
alg_set = [1:2];
var_set = [512,512+256,512+256*2,512+256*3,512+256*4];
M_set = [150,200];
Nsim = ones(1,length(var_set))*500;                                         % number of simulations
method = char('TS-OAMP-MA w/o sub-array', 'TS-OAMP-MA w/ sub-array');

AER = zeros(alg_type,length(var_set),length(M_set));                        % activity error rate
Pfa = zeros(alg_type,length(var_set),length(M_set));                        % false alarm probability
Pmd = zeros(alg_type,length(var_set),length(M_set));                        % miss detection probability
NMSE = zeros(alg_type,length(var_set),para.niter,length(M_set));            % NMSE

for i_M = 1:2
    para.M = M_set(i_M);
    for i_var=1:length(var_set)
        para.Nr = var_set(i_var);
        para.Num_sub = i_var+1;
        para.Nr_sub = para.Nr/para.Num_sub;

        %% Dictionnary Matrix
        sin_index = -1:2/para.Nr:1-2/para.Nr;
        D_far = exp(1i*pi*(0:para.Nr-1).' * sin_index)/sqrt(para.Nr);

        redunt = 2;
        sin_index = -1:2/redunt/para.Nr:1-2/redunt/para.Nr;
        D_rdt = exp(1i*pi*(0:para.Nr-1).' * sin_index)/sqrt(para.Nr*redunt);

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

        for i_sim = 1:Nsim(i_var)
            %% Channel Generation
            H = channel_generation_near_field_ULA_new(para);  

            %% Activity Generation
            Index_active = randperm(para.K,para.Ka);                         % actice device set
            act_flag = zeros(para.K,1);
            act_flag(Index_active) = 1;                                      % activity indicators
            Mask = zeros(para.K,para.Nr);
            for i_k = 1:para.Ka
                Mask(Index_active(i_k),:) = 1;
            end

            %% Measurement Matrix
            F = dftmtx(para.K)/sqrt(para.K);
            S = F(randperm(para.K,para.M),:);                                % pilot matrix

            %% Transmission
            H_space = H.*Mask;     
            Z= S * H_space;
            noise = sqrt(para.Pn/2) * (randn(size(Z)) + 1i*randn(size(Z)));
            Y = Z + noise;
            para.adj = sqrt(2)*0.675/3;                                      % to contorl the quantization range
            para.Ymax = para.adj *max(max(abs((Y)),[],'all'));
            Y_norm = Y / para.Ymax;                                          % normalization ( AGC)

            for alg_sel_ii = 1:length(alg_set)
                alg_sel = alg_set(alg_sel_ii);
                %% Quantization
                [Y_Q_R,~,~] = quantizer(real(Y_norm), para.Nbits, para.Nmax, -para.Nmax);                    % real part
                [Y_Q_I,delta_h, codebook_h] = quantizer(imag(Y_norm), para.Nbits, para.Nmax, -para.Nmax);    % iamginary part
                Y_Q = Y_Q_R+ 1i * Y_Q_I;                                     % complexy signal

                %% Add Mixed ADC
                para.Num_mix = para.Nr/ para.intvl;                          % total number of mixed-ADCs
                para.Num_mix_sub = para.Num_mix/para.Num_sub;
                mix_index = para.intvl/2:para.intvl:para.Nr;                 % uniformly distributed

                Y_wave = Y_Q;         
                Y_wave(:,mix_index) = Y_norm(:,mix_index);                   % high-resolution ADCs = infinite quantization

                %% Joint Activity Detection and Channel Estimation
                para.Nvar2 = para.Pn/para.Ymax/para.Ymax;                    % normalized noise variance
                switch alg_sel
                    case {1}  % TS-OAMP w/o sub-array
                        [h_est, act_hat, ~, v_s_h5, para] = TS_OAMP_withoutSubarray(Y_wave,S,codebook_h, para, mix_index, D_rdt);
                        h_est = h_est * para.Ymax;
                        H_est = h_est(:,:,para.niter);

                    case {2}  % TS-OAMP w/ sub-array
                        [h_est, act_hat, ~, v_s_h5, para] = TS_OAMP(Y_wave,S,codebook_h, para, mix_index, D_sub);   
                        h_est = h_est * para.Ymax;
                        H_est = h_est(:,:,para.niter);
                        
                    otherwise
                        disp('error! inexistent algorithm');
                end

                %% AER
                AER_tmp = sum(abs(act_hat - act_flag),'all')/para.K;
                Pfa_tmp = sum(abs(act_hat(find(act_flag==0)) - act_flag(find(act_flag==0))),'all')/(para.K-para.Ka);
                Pmd_tmp = sum(abs(act_hat(find(act_flag==1)) - act_flag(find(act_flag==1))),'all')/para.Ka;
                AER(alg_sel,i_var,i_M) = AER(alg_sel,i_var,i_M) + AER_tmp;
                Pfa(alg_sel,i_var,i_M) = Pfa(alg_sel,i_var,i_M) + Pfa_tmp;
                Pmd(alg_sel,i_var,i_M) = Pmd(alg_sel,i_var,i_M) + Pmd_tmp;

                %% NMSE
                for i_nmse=1:para.niter
                    h_est_tmp = h_est(:,:,i_nmse);
                    NMSE_temp = norm(h_est_tmp - reshape(H_space,para.K,[]),'fro')^2/norm(reshape(H_space,para.K,[]),'fro')^2;
                    NMSE(alg_sel, i_var, i_nmse,i_M) = NMSE(alg_sel, i_var, i_nmse,i_M) + NMSE_temp;
                end
                NMSE_temp = norm(H_est - reshape(H_space,para.K,[]),'fro')^2/norm(reshape(H_space,para.K,[]),'fro')^2;

                %% display
                if mod(i_sim,min(10,Nsim(i_var))) == 0 %|| i_sim == 1
                    fprintf('sim=%d, Pt=%d dBm, M=%d, Nr=%d, %s: AER=%.6f, AER_tmp=%.6f, NMSE=%.2f dB, NMSE_tmp=%.2f dB, Pfa=%.4f, Pmd=%.4f \n',...
                      i_sim, para.Ps_dBm, para.M, para.Nr, method(alg_sel,:),AER(alg_sel, i_var)/i_sim, AER_tmp,...
                      10*log10(NMSE(alg_sel, i_var, end)/i_sim),10*log10(NMSE_temp), Pfa(alg_sel, i_var)/i_sim, Pmd(alg_sel, i_var)/i_sim);
                end
            end

            if mod(i_sim, 200)==0
                Perform_tmp.NMSE = 10*log10(NMSE ./ repmat(i_sim,1,1,para.niter,1));
                Perform_tmp.AER = AER ./ i_sim;
                Perform_tmp.Pfa = Pfa ./ i_sim;
                Perform_tmp.Pmd = Pmd ./ i_sim;
                save Performance_tmp.mat Perform_tmp;
            end
        end
    end
end
Perform.NMSE = 10*log10(NMSE ./ repmat(Nsim,1,1,para.niter,length(M_set)));
Perform.AER = AER ./ Nsim;
Perform.Pfa = Pfa ./ Nsim;
Perform.Pmd = Pmd ./ Nsim;
disp('Finished all');
    
% save GainvsNr_Bit2.mat Perform;
%% Plot
var_set = [512,512+256,512+256*2,512+256*3,512+256*4];
MarkerSize = 6;
LineWidth = 2;
Fontsize = 15;
figure;
p1=plot(var_set,Perform.NMSE(1,:,end,1),'r-.s','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p2=plot(var_set,Perform.NMSE(2,:,end,1),'r-s','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p3=plot(var_set,Perform.NMSE(1,:,end,2),'b-.s','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
p4=plot(var_set,Perform.NMSE(2,:,end,2),'b-s','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);grid on;hold on;
set(gca,'FontSize',Fontsize);
legend('TS-OAMP w/o subarray, $M$ = 150','TS-OAMP, $M$ = 150','TS-OAMP w/o subarray, $M$ = 200','TS-OAMP, $M$ = 200','interpreter','latex');
xlabel('Number of antennas $N_r$','Fontsize',Fontsize,'interpreter','latex');
ylabel('NMSE','Fontsize',Fontsize,'interpreter','latex');
    


