function [lambdas,omegas,Sm,misclass,mse,time] = minimax_test_c_anis_blocks_exp_digits(sigma, block_sizes, mu, x, y, testx, testy, iters, ...
    tol_lambda, tol_opt, tol_e, tol_g, sigma0, lambda0, dc, dc_iters, dc_draws)

% sigma : #blocks x 2 matrix of [sigma_min, sigma_max] 
% sigma0 : n x #blocks matrix

[m,d] = size(x);

if (nargin == 14)
    dc = 0;
    dc_iters = 0;
end

tic;
[lambdas,omegas,K,Sm] = minimax_sq_c_exp_anis_blocks(sigma, block_sizes, mu, x, y, iters, tol_lambda, tol_opt, tol_e, tol_g, sigma0, ...
    lambda0, dc, dc_iters, dc_draws);
time = toc;

[mt,d] = size(testx);
B = length(block_sizes);
block_starts = [1,1+cumsum(block_sizes)];
for k=1:B
    Dxtx(k,:,:) = dist_matrix(x(:,block_starts(k):block_starts(k+1)-1),testx(:,block_starts(k):block_starts(k+1)-1));
end
Ktestx = zeros(m,mt);
for i=1:length(lambdas)
    Ktestx = Ktestx + lambdas(i)*exp(-squeeze(sum(Dxtx.*repmat(omegas(:,i),[1,m,mt]),1)));
end
    
f = y'*inv(K+mu*eye(m))*Ktestx;
misclass = length(find((sign(f')-testy)))/length(testy);
fprintf('minimax_test_c_anis_exp_digits() : \t %.8g error\n',misclass);
mse = sum((f'-testy).^2)/length(testy);