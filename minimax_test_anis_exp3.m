function [lambdas,Sm,err] = minimax_test_anis_exp3(x, y, testx, testy, mu, iters, tol_lambda, tol_g, recycle, lambdas0, ...
    sigma_l, sigma_u, n, blocks)

% Grid
% Learn with anisotropic Gaussian kernels K(x,y) = exp(-(x-y)'inv(cov)(x-y))
% cov: covariance matrices

[m,d] = size(x);
d1 = blocks(1);
d2 = blocks(2);
D1 = dist_matrix(x(:,1:d1),x(:,1:d1));
D2 = dist_matrix(x(:,d1+1:d1+d2),x(:,d1+1:d1+d2));

[lambdas,Sm] = minimax_sq2(D1, D2, mu, y, iters, tol_lambda, tol_g, recycle, lambdas0, sigma_l, sigma_u, n);
Kopt = zeros(m,m);
for i = 1:n
    for j = 1:n
        sigma1 = sigma_l+(i-1)*(sigma_u-sigma_l)/(n-1);
        sigma2 = sigma_l+(j-1)*(sigma_u-sigma_l)/(n-1);
        Kopt = Kopt + lambdas(n*(i-1)+j) * exp(-D1/(sigma1^2)) .* exp(-D2/(sigma2^2));
    end
end

[mt,d] = size(testx);
Dt1 = dist_matrix(x(:,1:d1),testx(:,1:d1));
Dt2 = dist_matrix(x(:,d1+1:d1+d2),testx(:,d1+1:d1+d2));

k = zeros(m,mt);
for i=1:n
    for j=1:n
        sigma1 = sigma_l+(i-1)*(sigma_u-sigma_l)/(n-1);
        sigma2 = sigma_l+(j-1)*(sigma_u-sigma_l)/(n-1);
        k = k + lambdas(n*(i-1)+j) * exp(-Dt1/(sigma1^2)) .* exp(-Dt2/(sigma2^2));
    end
end

f = kernel_regression(Kopt,y,mu)*k;
err = length(find(sign(f')-testy))/length(testy);
fprintf('Error = %f\n',err);
