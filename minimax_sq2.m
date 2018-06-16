function [lambdas,Sm] = minimax_sq2(D1, D2, mu, y, iters, tol_lambda, tol_g, recycle, lambdas0, sigma_l, sigma_u, n)

% Best convex combination of kernels for square loss 
% through iterative optimization wrt c and lambda
% K : n x m x m array of n kernel matrices

if (abs(sum(sum(lambdas0))-1) > eps)
    disp('Lambdas should sum to 1');
    return;
end
m = length(y);
lambdas = lambdas0;
Sm = inf;
j = 0;

for t=1:iters
disp(t);
    if ~recycle
        active = find(lambdas);
        if (length(active)==n^2)
            break;
        end
    end
    
    Kt = zeros(m,m);
    for i = 1:n
        for j = 1:n
            sigma1 = sigma_l+(i-1)*(sigma_u-sigma_l)/(n-1);
            sigma2 = sigma_l+(j-1)*(sigma_u-sigma_l)/(n-1);
            Kt = Kt + lambdas(n*(i-1)+j) * exp(-D1/(sigma1^2)) .* exp(-D2/(sigma2^2));
        end
    end
        
    Sm_prev = Sm;
    Sm = mu*y'*inv(Kt+ mu*eye(m))*y;
    if (abs(Sm-Sm_prev)/Sm < tol_g & j==n^2-1)
        break;
    end    
    
    ch = -inv(Kt+mu*eye(m))*y;
    g = ch'*Kt*ch;
    if recycle
        candidates = [1:n^2];
    else
        candidates = setdiff([1:n^2],active);
    end
    found = false;
    for j=candidates(randperm(length(candidates)))
        jj = rem(j,n);
        ii = (j-jj)/n+1;
        sigma1 = sigma_l+(ii-1)*(sigma_u-sigma_l)/(n-1);
        sigma2 = sigma_l+(jj-1)*(sigma_u-sigma_l)/(n-1);
        Kj = exp(-D1/(sigma1^2)) .* exp(-D2/(sigma2^2));
        if (ch'*Kj*ch > g)
            found = true;
            break;
        end
    end

    if (~found)         % Optimum reached (no kernel can improve)
        break;
    end
    
    % Search for lambda minimizing phi
    l_prev = inf;
    l = 0.5;
    while (abs(l-l_prev) > tol_lambda)
        tempi = inv(l*Kj+(1-l)*Kt+mu*eye(m));
        cKl = -2*mu*tempi*y;
        temp = (Kt-Kj)*cKl;
        dphi = cKl'*temp/(4*mu);
        ddphi = temp'*tempi*temp/(2*mu);
        l_prev = l;
        l = l-dphi/ddphi;
        if (l>=1)
            l=1;
            break;
        elseif (l<=0)
            l=0;
            break;
        end
    end
    
    % update kernel
    lambdas = l*[zeros(j-1,1);1;zeros(n^2-j,1)] + (1-l)*lambdas;
    
    %     % DEBUG
    %     figure;
    %     hold on;
    %     for l=0:0.005:1
    %         plot(l,-y'*inv(l*Kj+(1-l)*Kt+mu*eye(m))*y);
    %     end
    %     % END DEBUG
    templ(t,:) = lambdas';  
    tempS(t) = Sm;
end

%save('foo','templ','tempS');
disp(sprintf('Ended after %d iterations\n',t));
