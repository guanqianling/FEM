function z = FEM_1d_linear_Poisson(n)

if nargin < 1
    n = 16;
end

%% pre_setting 
format long e

solu = @u_real;
diff_u = @grad_u_real;
rhs = @f;
c = @coeff_c;
BC = [3, 1; 1, cos(1); 1, 0];

gs1_wt = [0.118463442528095, 0.239314335249683,...
          0.284444444444445, 0.239314335249683, 0.118463442528095];
gs1_pt = [0.046910077030668, 0.230765344947158, 0.5, ...
          0.769234655052842, 0.953089922969332];
gs1_n = 5;

A_pre = [1, -1; -1, 1] * n;

u_pre = [1 - gs1_pt; gs1_pt];

diff_pre = [-ones(1, gs1_n); ones(1, gs1_n)] * n;

%% stiffness_matrix
phy_gs = repmat(gs1_pt/n, n, 1) + repmat((0:1/n:1 - 1/n)', 1, gs1_n);

c_global = repmat(reshape(c(phy_gs)', 1, 1, gs1_n, n), 2, 2, 1, 1);

uv_global = repmat(A_pre, 1, 1, gs1_n, n);

ind = sort([1:n + 1 2:n]);
i = reshape(repmat(reshape(ind, 2, n), 2, 1), 2, 2 * n);
j = repmat(ind, 2, 1);
v = reshape(tensorprod(c_global .* uv_global, gs1_wt, 3, 2), 2, 2 * n);

A = sparse(i, j, v, n + 1, n + 1);

%% right_side
F_global = rhs(phy_gs);
F_global = F_global(reshape(ones(2, 1) * (1:n), 2 * n, 1), :);

v_global = repmat(u_pre, n, 1);

T = sparse(ind, 1:2 * n, ones(1, 2 * n), n + 1, 2 * n);
F = T * (F_global .* v_global) * gs1_wt' / n;

%% boundary_condition
if BC(1,1) == 1
    A(1, :) = [1, zeros(1, n)];
    F(1) = BC(2, 1);
else        
    F(1) = F(1) - c(0) * BC(3, 1);
    if BC(1, 1) == 3
        A(1, 1) = A(1, 1) - c(0) * BC(2, 1);
    end
end
if BC(1,2) == 1
    A(end, :) = [zeros(1, n), 1];
    F(end) = BC(2, 2);
else
    F(end) = F(end) + c(1) * BC(3, 2);
    if BC(1, 2) == 3
        A(end, end) = A(end, end) + c(1) * BC(2, 2);
    end
end

%% solve_AF
uh = A \ F;

%% error_estimate
u_real_global = solu(phy_gs);

T = sparse(sort([1:n 1:n]), 1:2 * n, ones(1, 2 * n), n ,2 * n);
u_global = T * (repmat(uh(ind), 1, gs1_n) .* repmat(u_pre, n, 1));

err = u_real_global - u_global;
err_L2 = sqrt(sum(err .* err * gs1_wt' / n));

err_Linf = max(max(abs(err)));

diff_u_real_global = diff_u(phy_gs);

diff_u_global = T * (repmat(uh(ind), 1, gs1_n) .* repmat(diff_pre, n, 1));

err = diff_u_real_global - diff_u_global;
err_H1_semi = sqrt(sum(err .* err * gs1_wt' / n));

z = [err_L2, err_Linf, err_H1_semi];
end

function z = coeff_c(x)
z = exp(x);
end

function z = f(x)
z = -exp(x) .* (cos(x) - 2 * sin(x) - x .* cos(x) - x .* sin(x));
end

function z = u_real(x)
z = x .* cos(x);
end

function z = grad_u_real(x)
z = cos(x) - x .* sin(x);
end