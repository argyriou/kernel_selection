function [x,fx] = cutting_plane_fmincon_vertex(g,h,a,y0,subg,subh,suba,S0,gmin,gmax,tol_opt,tol_lambda,maxit,draws,thres)

% Cutting plane method for DC programs of the form
% min f(x) = g(x)-h(x) s.t. a(x) <= 0 (see [Horst,Thoai])
% y0: a strictly feasible initial point (i.e. a(y0) < 0)
% subg: subgradient routine for g
% subh: subgradient routine for h
% suba: subgradient routine for a
% gmin: routine for minimizing convex function inside simplex
% gmax: routine for maximizing convex function inside simplex
% tol_opt : fractional tolerance for optimizations

% Hyperplane a'*x+b*t+c represented as [a;b;c]
% Polytope represented as
% a) (d+1) x #vertices matrix of vertices
% b) (d+2) x #constraints matrix of hyperplanes (constraints a'*x+b*t+c <= 0)
% c) index of active constraints per vertex -- #vertices x #constraints
% d) adjacency matrix for vertices
% S0: initial simplex containing feasible set X -- cell array of a,b,c,d
% (with dim = d-1)

    function ret = ht(x)
        ret = -feval(h,x(1:d))+x(d+1)';
    end
    function ret = f(x)
        ret = feval(g,x)-feval(h,x);
    end

if (nargin==14)
    thres = -inf;
end
thres

d = length(y0);
g_y0 = feval(g,y0);
w0 = g_y0-feval(h,y0);
wbar = feval(gmin,g,S0,tol_opt)-feval(gmax,h,S0,tol_opt);
tbar = feval(gmax,g,S0,tol_opt)-wbar+2*eps;
s = feval(subg,y0);

% Construct P0 using S0
cS0 = S0{2};
[foo,cnum] = size(cS0);
cP = [[cS0(1:d,:);zeros(1,cnum);cS0(d+1,:)],[zeros(d,1);1;-tbar],[s;-1;-s'*y0+g_y0-w0]];
VS0 = S0{1};
ht0 = min([-feval(h,VS0)+tbar ; -feval(h,VS0)+VS0'*s-repmat(y0',d+1,1)*s+g_y0-w0]);

% Bounds for polytope
vmin = [min(VS0,[],2);min(s'*VS0-s'*y0+g_y0-w0)];   % Bounds for polytope
vmax = [max(VS0,[],2);tbar];

yk = y0;
wk = w0;

k = 1;
pos_min = 0;
min_repeat = 0;
options = optimset('TolFun',tol_opt*abs(wk),'TolX',tol_opt*norm(y0),'TolCon',tol_opt*norm(y0),'Display','off',...
    'LargeScale','off','MaxFunEvals',100,'MaxIter',100);

while (k <= maxit)
    % Minimize -h(xk)+tk on the polytope vertices
    % First find random vertex using linprog

    % Find a better value by sampling
    A = cP(1:d+1,:)';
    C = A*vmin + cP(d+2,:)';
    [v0,simplex_flag] = call_simplex_init(A,C,eps,5*(2*cnum+2*k+d+5));
    if (~simplex_flag)
        continue;
    end
    v0 = v0 + vmin;

    [vtemp,ht_minval,exitflag,output,mults] = fmincon(@ht,v0,cP(1:d+1,:)',-cP(d+2,:)',[],[],[],[],[],options);

    if (exitflag < 1)
        fprintf('exitflag %d \n',exitflag);
        if (min_repeat < draws)
            min_repeat = min_repeat + 1;
        else
            x = yk;
            fx = wk;
            return;
        end
        continue;
    end
    min_repeat = 0;

    % vtemp might not be a vertex yet
    % Find a vertex on the face where vtemp lies
    active = find(mults.ineqlin);
    if (length(active)==d)
        A = [cP(1:d+1,:)';-cP(1:d+1,active)'];
        C = A*vmin+[cP(d+2,:)';-cP(d+2,active)'];
        [vk,simplex_flag] = call_simplex_init(A,C,eps,5*(2*cnum+2*k+3*d+5));
        if (~simplex_flag)
            disp('fmincon failed');
            continue;
        end
        vk = vk+vmin;
    else
        if (cP(1:d+1,:)'*vtemp+cP(d+2,:)' < eps)
            vk = vtemp;
        else
            disp('fmincon failed');
            continue;
        end
    end
    k = k+1;

    % If -h(x)+t >= 0 in the polytope then return

    if (ht_minval > tol_opt*ht0)
        if (pos_min == draws)
            x = yk;
            fx = wk;
            fprintf('cutting_plane() : \t Converged in %d iterations\n\n', k);
            return;
        else
            pos_min = pos_min + 1;
        end
    else
        pos_min = 0;
    end

    xk = vk(1:d);
    tk = vk(d+1);

    % Distinguish cases
    axk = feval(a,xk);
    if (axk <= 0)			% xk is feasible
        sk = feval(subg,xk);
        gxk = feval(g,xk);
        hxk = feval(h,xk);
        if (gxk-hxk < wk)
            [xk_min,fval,exitflag] = fmincon(@f,xk,cS0(1:d,:)',-cS0(d+1,:)',[],[],[],[],[],options);
            if (exitflag >=1 & gxk-hxk > fval)
                yk = xk_min;
                wk = fval;
            else
                yk = xk;
                wk = gxk-hxk;
            end
            if (wk < thres)
                x = yk;
                fx = wk;
                fprintf('cutting_plane() : \t Passed below threshold in %d iterations\n\n', k);
                return;
            end
        end
        lk = [sk;-1;gxk-wk-sk'*xk];	% New constraint
    else				% xk is not feasible
        beta2 = feval(g,xk)-tk-wk;
        if (axk > beta2)			% Get a subgradient of beta at (xk,tk)
            sk = [feval(suba,xk);0];
            beta_xktk = axk;
        else
            sk = [feval(subg,xk);-1];
            beta_xktk = beta2;
        end
        lk = [sk;beta_xktk-sk'*[xk;tk]];		% New constraint


        % Search for the zero of beta on the line segment (xk,tk),(y0,tbar) using Newton's method
        dl = -beta_xktk/(sk'*[xk-y0;tk-tbar]);
        lambda = 1+dl;
        newton_iters = 0;
        while (abs(dl) > tol_lambda & newton_iters < 50)
            newton_iters = newton_iters + 1;
            xl = lambda*xk+(1-lambda)*y0;
            tl = lambda*tk+(1-lambda)*tbar;
            beta1 = feval(a,xl);
            beta2 = feval(g,xl)-tl-wk;
            if (beta1 > beta2)
                sbeta_l = [feval(suba,xl);0];
                beta_l = beta1;
            else
                sbeta_l = [feval(subg,xl);-1];
                beta_l = beta2;
            end
            irate = sbeta_l'*[xk-y0;tk-tbar];
            dl = -beta_l/irate;
            if (abs(irate) < eps | abs(dl) < eps)
                break;
            end
            if (lambda+dl > 1 | lambda+dl < 0)
                break;
            end
            lambda = lambda + dl;
        end

        zk = lambda*xk+(1-lambda)*y0;
        gzk = feval(g,zk);
        hzk = feval(h,zk);
        if (gzk-hzk < wk)
            [zk_min,fval,exitflag] = fmincon(@f,zk,cS0(1:d,:)',-cS0(d+1,:)',[],[],[],[],[],options);
            if (exitflag >=1 & gzk-hzk > fval)
                yk = zk_min;
                wk = fval;
            else
                yk = zk;
                wk = gzk-hzk;
            end
            if (wk < thres)
                x = yk;
                fx = wk;
                fprintf('cutting_plane() : \t Passed below threshold in %d iterations\n\n', k);
                return;
            end
        end
    end

    % Update polytope by adding new constraint.
    cP = [cP,lk];
end

x = yk;
fx = wk;

fprintf('cutting_plane() : \t Could not converge. Reached within %.g of optimal value\n\n',ht_minval);

end
