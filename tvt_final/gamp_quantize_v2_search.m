function [xstore, Xhat_final, lambda_final, iter, iterIn, v, iter_final] = gamp_quantize_v2_search(Yq, Phi, rho, sigMax, sigMin, Nbits, niterIn, niterOut,...
                             nvar, damp, tol, tol_re, nns_sel, Nco, cell_num, d, H, para)
% AMP-based CS recovery with quantized measurements MMV
global y_res
global yq_res
global nmse_iter
y_res_min = +inf;
yq_res_min = +inf;
%% Initialization 
% problem dimension
[G, M, P] = size(Yq);
[~, N, ~] = size(Phi);
N_BS = M/cell_num;
Kc = N/cell_num;

% AMP for SLM
lambda = rho.*ones(N,M,P);
xmean = 0;
Xhat = xmean.*ones(N,M,P);
xvar = ones(N,M,P);
v = xvar;
V = ones(G,M,P);
Z = Yq;

% nearest neighboor sparsity pattern learning
index = ones(N,M,P);
index_left = [zeros(N,1,P), index(:,1:M-1,:)];
index_right = [index(:,2:M,:), zeros(N,1,P)];
index_ahead = cat(3, zeros(N,M,1), index(:,:,1:P-1));
index_latter = cat(3, index(:,:,2:P), zeros(N,M,1));
index_3D = index_left + index_right + index_ahead + index_latter; % index_up + index_down + 
clear index
clear index_left
clear index_right
clear index_ahead
clear index_latter

% Module B
lar_num = 1e6;
sma_num = 1e-6;
z_A_ext = zeros(G,M);
% v_A_ext = lar_num;
v_A_ext = trace(Phi*Phi')/G*rho;
delta = (sigMax-sigMin)/2^Nbits;

% allocate memory
C = zeros(N,M,P);
D = zeros(N,M,P);
L_cal = zeros(N,M,P);
pai = zeros(N,M,P);
A = zeros(N,M,P);
B = zeros(N,M,P);
Yhat = zeros(G,M,P);
xstore = zeros(N,M,niterIn*niterOut);
iter_iii=1;

%% Turbo iterations
NMSE_yq = zeros(niterOut,1);
for iter = 1:niterOut
    
    Xhat_pre_out = Xhat;
    lambda_pre_out = lambda;
    
    % Module B: component-wise MMSE estimate of z = Ax
    z_B_post = zeros(G,cell_num*N_BS);
    v_B_post = 0;
    for c = 1:cell_num
        idAnt = (c-1)*N_BS+1:c*N_BS;
        delta = (sigMax(c)-sigMin(c))/2^Nbits;
        [zr_tmp, vr_tmp] = mmse_estimator(real(Yq(:,idAnt)), real(z_A_ext(:,idAnt)), v_A_ext.*ones(G,N_BS),...
                                              nvar, Nbits, delta, para);
        [zi_tmp, vi_tmp] = mmse_estimator(imag(Yq(:,idAnt)), imag(z_A_ext(:,idAnt)), v_A_ext.*ones(G,N_BS),...
                                              nvar, Nbits, delta, para);
        z_B_post(:,idAnt) = zr_tmp + 1i.* zi_tmp;
        v_B_post = v_B_post + mean(mean(vr_tmp)) + mean(mean(vi_tmp)); 
    end
    v_B_post = v_B_post /cell_num;
%     [z_B_post_r, v_B_post_r] = mmse_estimator(real(Yq), real(z_A_ext), v_A_ext.*ones(G,M),...
%                                               nvar, Nbits, delta);
%     [z_B_post_i, v_B_post_i] = mmse_estimator(imag(Yq), imag(z_A_ext), v_A_ext.*ones(G,M),...
%                                               nvar, Nbits, delta);
%     z_B_post = z_B_post_r + 1i.* z_B_post_i;
%     v_B_post = mean(mean(v_B_post_r)) + mean(mean(v_B_post_i)); 
    
    % extrinct imformation of module B
    sigma2_tilde = v_B_post.*v_A_ext ./ (v_A_ext-v_B_post);   
    sigma2_tilde = lar_num*(sigma2_tilde<0) + sigma2_tilde.*(sigma2_tilde>0);
    sigma2_tilde = min(sigma2_tilde,lar_num);
    sigma2_tilde = max(sigma2_tilde,sma_num);
    y_tilde = sigma2_tilde.*(z_B_post./v_B_post - z_A_ext./v_A_ext);    
    
    %%  Module A: AMP iteration for SLM 
%     [Xhat, lambda, iter, v] = gmmv_amp(y_tilde, Phi, damp, niterIn,...
%                          tol, tol, nns_sel, Nco, cell_num, rho, varH, nvar, d, H_re); 
%     [Xhat, lambda, ~, v] = gmmv_amp(y_tilde, Phi, damp, niterIn,...
%                          tol, tol, 0, 1, 1, 0.05, ones(N,M), nvar, ones(1,N), ones(1,N)); 
    
    % inner iteration
    NMSE_yq_sum = 0;
    for iterIn = 1:niterIn
        
        % initialization
%         lambda = rho.*ones(N,M,P);
%         Z = y_tilde;

        Xhat_pre = Xhat;
        V_pre = V;
        
        for p = 1:P
            % factor node update
            V(:,:,p) = damp.*V_pre(:,:,p) + (1-damp).*abs(Phi(:,:,p)).^2*v(:,:,p);
            Z(:,:,p) = damp.*Z(:,:,p) + (1-damp).*(Phi(:,:,p)*Xhat(:,:,p) - ...
                                  V(:,:,p).*(y_tilde(:,:,p)-Z(:,:,p))./(sigma2_tilde+V_pre(:,:,p)));
     
            % variable node update
            D(:,:,p) = 1 ./ ((abs(Phi(:,:,p)).^2).' * (1./(sigma2_tilde+V(:,:,p))));
            C(:,:,p) = Xhat(:,:,p) + D(:,:,p).*(Phi(:,:,p)'*((y_tilde(:,:,p)-Z(:,:,p))./(sigma2_tilde+V(:,:,p))));

            % posterior mean and variance
            L_cal(:,:,p) =0.5.*(log(D(:,:,p)./(D(:,:,p)+xvar(:,:,p))) + abs(C(:,:,p)).^2./D(:,:,p)...
                                                   - abs(C(:,:,p)-xmean).^2./(D(:,:,p)+xvar(:,:,p))); 
            pai(:,:,p) = lambda(:,:,p) ./ (lambda(:,:,p) + (1-lambda(:,:,p)).*exp(-L_cal(:,:,p))); 
            A(:,:,p) = (xvar(:,:,p).*C(:,:,p)+xmean.*D(:,:,p)) ./ (D(:,:,p)+xvar(:,:,p));
            B(:,:,p) = (xvar(:,:,p).*D(:,:,p)) ./ (xvar(:,:,p)+D(:,:,p));
            Xhat(:,:,p) = pai(:,:,p).*A(:,:,p);
            v(:,:,p) = pai(:,:,p).*(abs(A(:,:,p)).^2 + B(:,:,p)) - abs(Xhat(:,:,p)).^2;

            % reconstruct received signal
            Yhat(:,:,p) = Phi(:,:,p)*Xhat(:,:,p);
        end
        xstore(:,:,iter*niterIn+iterIn-niterIn)=Xhat;
%         % EM-based parameter learning
%         xmean = repmat(sum(pai.*A,1)./sum(pai,1), N, 1);
%         xvar = pai.*(abs(xmean-A).^2+B)./pai; 
        xvar = ones(size(pai)).*sum(pai.*(abs(xmean-A).^2+B),'all')./sum(pai,'all'); 
% if nns_sel ==4
        nvar_temp = abs(y_tilde-Z).^2./abs(1+V./sigma2_tilde).^2 + V./(1+V./sigma2_tilde);
        sigma2_tilde = sum(nvar_temp(:))/G/M/P;
%         nvar = sum(nvar_temp(:))/G/M/P;
% end

   
       % exploiting the structured sparsity for refining the update rule of sparsity ratio
       if nns_sel == 0  % each AP performs AUD and CE locally
          for nn = 1:cell_num
              pai_avg = sum(sum(pai((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :),2),3) ./ (N_BS*P);
              lambda((nn-1)*Kc+1:nn*Kc, (nn-1)*N_BS+1:nn*N_BS, :) = pai_avg(:, ones(1,N_BS), ones(1,P));
          end

       elseif nns_sel == 1  % Cloud computing  
             bi_ap = zeros(N,cell_num); % average at AP
             for c = 1:cell_num                                         
                 bi_ap(:,c) = sum(sum(pai(:, (c-1)*N_BS+1:c*N_BS, :),2),3) ./ (N_BS*P);
             end

             for nn = 1:N
                 d_sum = sum(1./d(:,nn),1);
                 bi_tmp = 0;
                 for c = 1:cell_num
                     bi_tmp = bi_tmp + 1/(d(c,nn)*d_sum)*bi_ap(nn,c);
                 end
                 lambda(nn, :, :) = bi_tmp;
             end

       elseif nns_sel == 2  %  Edge computing new
             bi_ap = zeros(N,cell_num); % average for each APs
             for c = 1:cell_num                                         
                 bi_ap(:,c) = sum(sum(pai(:, (c-1)*N_BS+1:c*N_BS, :),2),3) ./ (N_BS*P);
             end

             for nn = 1:N
                 [~,sel_cell] = sort(d(:,nn));
                 sel_cell = sel_cell(1:Nco);
                 d_sum = sum(1./d(sel_cell,nn),1);
                 bi_tmp = 0;
                 for c = 1:numel(sel_cell)
                     bi_tmp = bi_tmp + 1/(d(sel_cell(c),nn)*d_sum)*bi_ap(nn,sel_cell(c));
                 end
                  lambda(nn, :, :) = bi_tmp;
             end

       elseif nns_sel == 3  %  Edge computing
          for nn = 1:N
              [~,sel_cell] = sort(d(:,nn));
              sel_cell = sel_cell(1:Nco);

              pai_temp = 0;
              for c = 1:length(sel_cell)
                  pai_avg = sum(pai(nn,(sel_cell(c)-1)*N_BS+1:sel_cell(c)*N_BS))/N_BS;
                  pai_temp = pai_temp + pai_avg;
              end

              for c = 1:length(sel_cell)
                  lambda(nn,(sel_cell(c)-1)*N_BS+1:sel_cell(c)*N_BS) = pai_temp/length(sel_cell);
              end
           end

       else  % clustered sparsity 
%           pai_left = [zeros(N,1,P), pai(:,1:M-1,:)];
%           pai_right = [pai(:,2:M,:), zeros(N,1,P)];
%           pai_ahead = cat(3, zeros(N,M,1), pai(:,:,1:P-1));
%           pai_latter = cat(3, pai(:,:,2:P), zeros(N,M,1));
%           lambda = (pai_left + pai_right + pai_ahead + pai_latter) ./ index_3D;  
          
         % Common Sparsity     
         lambda = mean(pai,2);
         lambda = repmat(lambda, 1, M);
         
%          lambda = pai;                          % Sparsity        
%          cube=ones(N,M);
%          t0=zeros(N,M,2);
%          for i_nnspl=1:2
%              t1=cat(2,cube(:,2:end),zeros(N,1));               
%              t2=cat(2,zeros(N,1),cube(:,1:end-1));
% %              t3=cat(2,cube(:,3:end),zeros(K,2));               
% %              t4=cat(2,zeros(K,2),cube(:,1:end-2));
% %              t5=cat(2,cube(:,4:end),zeros(K,3));               
% %              t6=cat(2,zeros(K,3),cube(:,1:end-3));
%              t0(:,:,i_nnspl)=t1+t2;%+t3+t4+t5+t6;
%              cube=lambda;
%          end
%          lambda = t0(:,:,2)./t0(:,:,1); 
       end

       % stopping criteria
       NMSE_iter =  norm(Xhat(:)-Xhat_pre(:))^2 / norm(Xhat_pre(:))^2;
       NMSE_re = norm(y_tilde(:)-Yhat(:))^2 / norm(y_tilde(:))^2;
       NMSE_qre = norm(Yq(:)-Yhat(:))^2 / norm(Yq(:))^2;
       if NMSE_iter < tol || NMSE_re < tol_re
           break;           
       end
       
       NMSE_yq_sum = NMSE_yq_sum + NMSE_qre;
   
       % print status
       NMSE = norm(Xhat(:)-H(:))^2 / norm(H(:))^2;
       %y_res = norm(Xhat(:)-H(:))^2 / norm(H(:))^2;
%        fprintf('iterOut = %d, iterIn = %d, NMSE = %7.9f, y_res = %7.9f, yq_res = %7.9f, sigma2_tilde = %7.9f, v_A_ext = %7.9f,  v_B_post = %7.9f\n',...
%                 iter, iterIn, NMSE, NMSE_re,NMSE_qre,sigma2_tilde, v_A_ext, v_B_post);
    
            
       y_res(iter_iii) = NMSE_re;
       nmse_iter(iter_iii) = NMSE;
            
       if y_res_min > y_res(iter_iii)
          y_res_min = y_res(iter_iii);
          Xhat_final = Xhat;
%           lambda_final = lambda;
       end
       iter_iii=iter_iii+1;
    end
    
    yq_res(iter) = NMSE_yq_sum/iterIn;
    if yq_res_min > yq_res(iter)
       yq_res_min = yq_res(iter);
       Xhat_final = Xhat;
       lambda_final = pai;%lambda;
       iter_final = (iter-1)*niterIn+iterIn;
    else
       Xhat_final = Xhat_pre_out;
       lambda_final = lambda_pre_out;
       iter_final = (iter-1)*niterIn+iterIn;
       break;
    end
    
    % extrinct imformation of module A
    ext_mode = 1;
    if ext_mode == 0
        z_A_post = Phi*Xhat;
        v_A_post = mean(mean(abs(Phi).^2*v));
        v_A_ext = v_A_post.*sigma2_tilde ./ (sigma2_tilde-v_A_post);
    else
        v_A_ext = mean(mean(abs(Phi).^2*v));
    end
    v_A_ext = lar_num.*(v_A_ext<0) + v_A_ext.*(v_A_ext>0);
    v_A_ext = min(v_A_ext,lar_num);
    v_A_ext = max(v_A_ext,sma_num);
    if ext_mode == 0
        z_A_ext = v_A_ext.*(z_A_post./v_A_post - y_tilde./sigma2_tilde);    
    else
        z_A_ext = Phi*Xhat - V.*(y_tilde-Z)./(sigma2_tilde+V_pre);
    end
    
    NMSE_iter_out =  norm(Xhat(:)-Xhat_pre_out(:))^2 / norm(Xhat_pre_out(:))^2;
    NMSE_re = norm(y_tilde(:)-Yhat(:))^2 / norm(y_tilde(:))^2;
%     if NMSE_iter_out < tol || NMSE_re < tol_re 
%        break;
%     end



end
end