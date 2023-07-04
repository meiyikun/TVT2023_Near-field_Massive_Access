function [miu_storage, lambda, thegma2, v_s] = OAMP_MMV(Y_wave, S, codebook, para, D_far)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_wave: Received signals with mixed-resolution quantization
% S: Pilot matrix
% codebook : Quantization codebook
% para: System parameters
% D_sub: Traditional dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
M = size(Y_wave,1);     Nr = size(Y_wave,2);    K = size(S,2);
damping = para.damping;
niter = para.niter;
delta = codebook(2) - codebook(1);
v_delta = 1e-6;

% EM initialization
c = -10:0.01:10;                       c = exp(c);
Phinorm = 1-0.5*erfc(-c/sqrt(2));      phinorm = 1/sqrt(2*pi)*exp(-c.^2/2);
templambda0 = (1+c.^2).*Phinorm-c.*phinorm;
labmda0 = M/K*max((1-2*K*(templambda0)/M)./(1+c.^2-2.*(templambda0)));
lambda = repmat(labmda0,K,Nr);
thegma2 = norm(Y_wave,2)^2/101/M/Nr;
% thegma2 = para.Nvar2;
Psi=zeros(1,Nr);
for j=1:Nr
    Psi(j)=(norm(Y_wave(:,j),2)^2-M*thegma2)/norm(S,'fro')^2/labmda0;
end

v = ones(Nr,1);         u = zeros(K,Nr);
v_prev = zeros(Nr,1);   u_prev = zeros(K,Nr);
% vAext_prev = zeros(Nr,1);   zAext_prev = zeros(M,Nr);
nu_wave = ones(Nr,1);          u_wave = zeros(M,Nr);
v_s = zeros(niter, 1);
miu_storage = zeros(K,Nr,niter);
y = [reshape(real(Y_wave),[],1),reshape(imag(Y_wave),[],1)];
nu_ = zeros(Nr,2);     z_hat_ = zeros(M,Nr,2);
u_ = zeros(M,Nr,2);
for iter=1:niter
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
        yita1 = (sy.* u_(:,:,rti) - min(rb_1,rb)) ./ sqrt((thegma2 + nu_wave.')/2);
        yita2 = (sy.* u_(:,:,rti) - max(rb_1,rb)) ./ sqrt((thegma2 + nu_wave.')/2);
        
        temp1 = normpdf(yita1) - normpdf(yita2);
        temp2 = normcdf(yita1) - normcdf(yita2);
        temp3 = yita1.*normpdf(yita1) - yita2.*normpdf(yita2);
        temp4 = temp1./temp2;
        temp5 = temp3./temp2 + (temp1./temp2).^2;
        Small_toc = 1e-50;
        temp4(abs(temp2)<Small_toc) = - yita1(abs(temp2)<Small_toc);
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
    
    %% SLM
    % Transformation to angular-domain
    z_ext = z_ext*D_far;
    nu_ext = D_far'*diag(nu_ext)*D_far;
    nu_ext = real(diag(nu_ext));
    
    v=zeros(Nr,1);
    for j = 1:Nr
        v(j) = norm(z_ext(:,j) - S*u(:,j) ,'fro')^2/M - nu_ext(j);
    end
    v = max(v, 1e-20);
    
    r = u+K/M*(S'*(z_ext-S*u));
    tao = (K-M)/M*v + K/M*nu_ext;
    
    meanX = (Psi./(tao.'+Psi)).*r;
    varX = tao.'.*Psi./(tao.'+Psi);
    exponent = -(Psi./(Psi+tao.')./tao.').*abs(r).^2;
    
    pai_temp = max((1-lambda)./lambda.*((tao.'+Psi)./tao.').*exp(exponent), 1e-15);
    pai = 1./(1+pai_temp);
    miu = pai.*meanX;
    miu_storage(:,:,iter) = miu;
    gamma = pai.*varX+pai.*(1-pai).*abs(meanX).^2;
    
    gamma_mean = (mean(gamma,1)).';
    v = 1./(1./gamma_mean-1./tao);
    v = min(v,1e5);
    v = max(v,v_delta);
    u = (v.').*((miu./(gamma_mean.')-r./(tao.')));
    
    u = (1-damping)*u+damping*u_prev;   v = (1-damping)*v+damping*v_prev;
    u_prev = u;                         v_prev = v;
    
    %% EM update
    Psi = sum(pai.*(varX+abs(meanX).^2))./sum(pai);                         % Gaussian Variance of Signal
    
    lambda = reshape(pai, K, Nr);                                           % Sparsity
    cube = ones(K,Nr);
    t0 = zeros(K,Nr,2);
    for i_nnspl = 1:2
        t1 = cat(2,cube(:,2:end),zeros(K,1));
        t2 = cat(2,zeros(K,1),cube(:,1:end-1));
        t0(:,:,i_nnspl) = t1+t2;
        cube = lambda;
    end
    lambda = t0(:,:,2)./t0(:,:,1);
    
    thegma2 = mean(abs(z_ext-S*miu).^2+abs(S).^2*gamma,'all').';            % Noise variance
    
    u_wave = S*u*D_far';
    nu_wave = D_far*diag(v)*D_far';
    nu_wave = real(diag(nu_wave));
end
end



