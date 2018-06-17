load 'data/digitsoe';
load 'data/digits47';
load 'data/digits38';

n = 5;
d1 = 392;
d2 = 392;
mu = 1e-7;
sigma_u = 5000;
sigma_l = 500;

lambdas0 = ones(n^2,1)/(n^2);
lambdas0(n^2) = 1-sum(lambdas0(1:n^2-1));

[lambdas, Sm, misclassoe] = minimax_test_anis_exp3(xoe(1:500,:), yoe(1:500), txoe, tyoe, mu, 10, 1e-6, 1e-6, 1, lambdas0, sigma_l, sigma_u, n, [d1,d2]);
[lambdas, Sm, misclass47] = minimax_test_anis_exp3(x47(1:500,:), y47(1:500), tx47, ty47, mu, 100, 1e-6, 1e-6, 1, lambdas0, sigma_l, sigma_u, n, [d1,d2]);
[lambdas, Sm, misclass38] = minimax_test_anis_exp3(x38(1:500,:), y38(1:500), tx38, ty38, mu, 100, 1e-6, 1e-6, 1, lambdas0, sigma_l, sigma_u, n, [d1,d2]);

save(sprintf('results_grid_s%d_%d_%dx%d',sigma_l,sigma_u,n,n), 'misclassoe', 'misclass38','misclass47'); 

