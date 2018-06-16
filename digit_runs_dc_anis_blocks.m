load '~/data/colt/digitsoe';
load '~/data/colt/digits47';
load '~/data/colt/digits38';

n = 5;
d = 784;
m = 500;

mu = 1e-7;

sizes = {[784],[392,392],[261,261,262],[196,196,196,196]};

for i=1:4

    block_sizes = sizes{i};

    ranges = {[75,25000],[100,1e4],[500,5000]};
    
    for j=1:3

        range = ranges{j};
        sigma_u = range(2);
        sigma_l = range(1);
        sigma0 = repmat(sigma_l:(sigma_u-sigma_l)/(n-1):sigma_u,length(block_sizes),1)';
        lambda0 = ones(1,n)/n;

        sigma = repmat([sigma_u,sigma_l],length(block_sizes),1);
        maxit = 100;
        dcmaxit = 300;
        dc_draws = 10;
        tol_lambda = 1e-6;
        tol_opt = 1e-3; 
        tol_Q = 1e-12;
        tol_fun = 1e-6; 

        [lambdas,omegas,Sm,misclass47,mse,time47] = minimax_test_c_anis_blocks_exp_digits(sigma, block_sizes, mu, x47(1:m,:), y47(1:m), ...
            tx47, ty47, maxit, tol_lambda, tol_opt, tol_Q, tol_fun, sigma0, lambda0, 1, dcmaxit, dc_draws);
        save(sprintf('results_47_b%d_%f_%f.mat',i,sigma_l,sigma_u),'lambdas','omegas','Sm','misclass47','time47','sigma_u','sigma_l','block_sizes');

        [lambdas,omegas,Sm,misclassoe,mse,timeoe] = minimax_test_c_anis_blocks_exp_digits(sigma, block_sizes, mu, xoe(1:m,:), yoe(1:m), ...
            txoe, tyoe, maxit, tol_lambda, tol_opt, tol_Q, tol_fun, sigma0, lambda0, 1, dcmaxit, dc_draws);
        save(sprintf('results_oe_b%d_%f_%f.mat',i,sigma_l,sigma_u),'lambdas','omegas','Sm','misclassoe','timeoe','sigma_u','sigma_l','block_sizes');

        [lambdas,omegas,Sm,misclass38,mse,time38] = minimax_test_c_anis_blocks_exp_digits(sigma, block_sizes, mu, x38(1:m,:), y38(1:m), ...
            tx38, ty38, maxit, tol_lambda, tol_opt, tol_Q, tol_fun, sigma0, lambda0, 1, dcmaxit, dc_draws);
        save(sprintf('results_38_b%d_%f_%f.mat',i,sigma_l,sigma_u),'lambdas','omegas','Sm','misclass38','time38','sigma_u','sigma_l','block_sizes');

    end
end


