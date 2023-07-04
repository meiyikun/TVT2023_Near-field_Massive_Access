function [miu_storage, act_hat, thegma2, v_s, para] = TS_OAMP_withoutSubarray(Y_wave, S, codebook, para, mix_index, D_rdt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_wave: Received signals with mixed-resolution quantization
% S: Pilot matrix
% codebook : Quantization codebook
% para: System parameters
% mix_index: Index set of antennas with high-resolution ADCs
% D_red: Traditional dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization 
M = size(Y_wave,1);     Nr = size(Y_wave,2);    K = size(S,2);
damping = para.damping;
niter = para.niter;
Num_sub = 1;
Nr_sub = Nr/Num_sub;
Gr = size(D_rdt,2);
beta = Gr/Nr;
delta = codebook(2) - codebook(1); 
v_delta = 1e-6; 

%EM initialization
c = -10:0.01:10;                       c = exp(c);
Phinorm = 1-0.5*erfc(-c/sqrt(2));      phinorm = 1/sqrt(2*pi)*exp(-c.^2/2);
templambda0 = (1+c.^2).*Phinorm-c.*phinorm;
labmda0 = M/K*max((1-2*K*(templambda0)/M)./(1+c.^2-2.*(templambda0)));
lambda_s = repmat(labmda0,K,Nr);
lambda_a = repmat(labmda0,K,Gr);
thegma2 = norm(Y_wave(:,mix_index),'fro')^2/101/M/length(mix_index);
% thegma2 = para.Nvar2;
Psi = zeros(1,Nr);
for j=1:Nr
    Psi(j) = (norm(Y_wave(:,j),2)^2-M*thegma2)/norm(S,'fro')^2/labmda0;
end

v = ones(Nr,1);            u = zeros(K,Nr);
v_prev_s  = zeros(Nr,1);   u_prev_s  = zeros(K,Nr); 
v_prev_a = zeros(Nr,1);    u_prev_a = zeros(K,Nr); 
nu_wave = ones(Nr,1);      u_wave = zeros(M,Nr);
v_s = zeros(niter, 1);
miu_storage = zeros(K,Nr,niter);      
y = [reshape(real(Y_wave),[],1),reshape(imag(Y_wave),[],1)]; 
Y_wave_0 = reshape(Y_wave, M, Nr);
nu_ = zeros(Nr,2);     z_hat_ = zeros(M,Nr,2);
u_ = zeros(M,Nr,2);      
for iter=1:niter
    if iter<=para.T1
        %% Stage 1
        %% GLM, de-quantization
        u_(:,:,1) = real(u_wave);           u_(:,:,2) = imag(u_wave);
        for rti = 1:2
            y_ = y(:, rti);
            sy = sign(reshape(y_, M, Nr));
            r_low = y_ - delta/2;      r_up = y_ + delta/2;
            r_low(r_low <= codebook(1)) = -1/para.adj/sqrt(2);
            r_up(r_up >= codebook(end)) = 1/para.adj/sqrt(2);
            rb_1 = reshape(abs(r_low),M,Nr);
            rb = reshape(abs(r_up),M,Nr);
            yita = (sy.* u_(:,:,rti) - min(rb_1,rb)) ./ sqrt((thegma2 + nu_wave.')/2);
            xi = (sy.* u_(:,:,rti) - max(rb_1,rb)) ./ sqrt((thegma2 + nu_wave.')/2);
            
            temp1 = normpdf(yita) - normpdf(xi);
            temp2 = normcdf(yita) - normcdf(xi);
            temp3 = yita.*normpdf(yita) - xi.*normpdf(xi);
            temp4 = temp1./temp2;
            temp5 = temp3./temp2 + (temp1./temp2).^2;
            Small_toc = 1e-50;
            temp4(abs(temp2)<Small_toc) = - yita(abs(temp2)<Small_toc);
            temp5(abs(temp2)<Small_toc) = 1;
            
            z_hat_(:,:,rti) = u_(:,:,rti) + sy .* nu_wave.' ./ (sqrt(2*(thegma2 + nu_wave.'))).* temp4;
            nu_(:,rti) = mean(nu_wave.' / 2 - nu_wave.'.^2./(2*(thegma2+nu_wave.')).* temp5, 1).';
        end
        
        % extrinsic information
        nu = nu_(:,1) + nu_(:,2);
        z_hat = z_hat_(:,:,1) + 1i.*z_hat_(:,:,2);
        
        nu_ext = 1./(1./nu - 1./nu_wave);
        nu_ext = min(nu_ext,1e5);           nu_ext = max(nu_ext,v_delta);
        z_ext = nu_ext.'.* ((z_hat./nu.' - u_wave./nu_wave.'));
        
        % mixed-ADC
        z_ext_mixed = reshape(z_ext, M, Nr);
        z_ext_mixed(:,mix_index) = Y_wave_0(:,mix_index);
        z_ext = reshape(z_ext_mixed, M, Nr);
        
        nu_ext_mixed = nu_ext;
        nu_ext_mixed(mix_index) = thegma2;
        nu_ext = reshape(nu_ext_mixed, [], 1);
        
        %% SLM, spatial,domain
        v = zeros(Nr,1);
        for j = 1:Nr
            v(j) = norm(z_ext(:,j) - S*u(:,j) ,'fro')^2/M - nu_ext(j);
        end
        v = max(v, 1e-20);
        
        r = u+K/M*(S'*(z_ext-S*u));  
        tao = (K-M)/M*v + K/M*nu_ext;      
        
        meanX = (Psi./(tao.'+Psi)).*r;
        varX = tao.'.*Psi./(tao.'+Psi);
        exponent = -(Psi./(Psi+tao.')./tao.').*abs(r).^2;
        
        pai_temp = max((1-lambda_s)./lambda_s.*((tao.'+Psi)./tao.').*exp(exponent), 1e-15);
        pai = 1./(1+pai_temp);
        miu = pai.*meanX;
        miu_storage(:,:,iter) = miu;
        gamma = pai.*varX+pai.*(1-pai).*abs(meanX).^2;
        
        gamma_mean = (mean(gamma,1)).'; 
%         gamma_mean = (mean(gamma,'all')).'; 
        v = 1./(1./gamma_mean-1./tao);
        v = min(v,1e5);
        v = max(v,v_delta);
        u = (v.').*((miu./(gamma_mean.')-r./(tao.')));
        
        u = (1-damping)*u+damping*u_prev_s;   v = (1-damping)*v+damping*v_prev_s;
        u_prev_s = u;                         v_prev_s = v;
        
        % EM update
        Psi=sum(pai.*(varX+abs(meanX).^2),1)./sum(pai,1);                   % Gaussian Variance of Signal
        
        lambda_s = pai;                                                     % Sparsity
        if  thegma2 >= norm(Y_wave(:,mix_index),'fro')^2/11/M/length(mix_index)
            lambda_s = repmat(mean(lambda_s,2),1,Nr);
        else
            lambda_s = repmat(mean(lambda_s(:,mix_index),2),1,Nr);
        end
        
        thegma2 = mean(abs(Y_wave(:,mix_index)-S*miu(:,mix_index)).^2+gamma_mean(mix_index).','all').'; % Noise variance
    else
        %% Stage 2
        % Initialization 
        rou_final = reshape(lambda_s,K,[]);
        act_hat = (mean(rou_final,2)>para.threshold);
        if iter == (para.T1+1)
            u = u*D_rdt;
            u_prev_a = u_prev_s*D_rdt;
            v_prev_a = D_rdt'*diag(v_prev_s)*D_rdt;
            v_prev_a = real(diag(v_prev_a));

            Psi = D_rdt'*diag(Psi)*D_rdt;
            Psi = real(diag(Psi))';

            z_ext = z_ext*D_rdt; 
            nu_ext = D_rdt'*diag(nu_ext)*D_rdt;
            nu_ext = real(diag(nu_ext)); 
        end
        
        %% SLM, angular domain
        v = zeros(Gr,1);
        for j = 1:Gr
            v(j) = norm(z_ext(:,j) - S*u(:,j) ,'fro')^2/M - nu_ext(j);
        end
        v = mean(v);
        
        r = u+K/M*(S'*(z_ext-S*u));  
        tao = (K-M)/M*v + K/M*nu_ext;     
 
        meanX = (Psi./(tao.'+Psi)).*r;
        varX = tao.'.*Psi./(tao.'+Psi);
        exponent = -(Psi./(Psi+tao.')./tao.').*abs(r).^2;
        lambda_a = lambda_a .* act_hat;
        pai_temp = max((1-lambda_a)./lambda_a.*((tao.'+Psi)./tao.').*exp(exponent), 1e-15);
        pai = 1./(1+pai_temp);
        miu = pai.*meanX;
        miu_storage(:,:,iter) = miu*D_rdt';                                 
        gamma = pai.*varX+pai.*(1-pai).*abs(meanX).^2;   

         gamma_mean = (mean(gamma,1)).';  
%          gamma_mean = (mean(gamma,'all')).'; 
         v = 1./(1./gamma_mean-1./tao);
         v = min(v,1e5);  
         v = max(v,v_delta);
         u = (v.').*((miu./(gamma_mean.')-r./(tao.')));

         u = (1-damping)*u+damping*u_prev_a;   v = (1-damping)*v+damping*v_prev_a;
         u_prev_a = u;                         v_prev_a = v;  

         % EM update
         Psi = sum(pai.*(varX+abs(meanX).^2),1)./sum(pai,1);                % Gaussian Variance of Signal
         
         lambda_a = reshape(pai, K, beta*Nr_sub, Num_sub);                  % Sparsity    
         cube = ones(K,beta*Nr_sub, Num_sub);
         t0 = zeros(K,beta*Nr_sub, Num_sub,beta);
         for i_nnspl = 1:2
             t1 = cat(2,cube(:,2:end,:),zeros(K,1,Num_sub));               
             t2 = cat(2,zeros(K,1,Num_sub),cube(:,1:end-1,:));
             t0(:,:,:,i_nnspl) = t1+t2;
            cube = lambda_a;
        end
        lambda_a = t0(:,:,:,2)./t0(:,:,:,1); 
        lambda_a = reshape(lambda_a, K, Gr);
    end
        
    if iter <= para.T1
        u_wave = S*u; 
        nu_wave = v;
    end
end
end



