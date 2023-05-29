function boundary_condition_R(T)

global p s
global C F
global dof 
global coeff_c bd_r bd_q
global gs1_pt gs1_wt gs1_n

bound_idx = find(s(4, :) == -3);
n = length(bound_idx);

u_global = reshape([1 - gs1_pt; gs1_pt], 2, 1, gs1_n);
v_global = reshape([1 - gs1_pt; gs1_pt], 1, 2, gs1_n);

uv_global = repmat(pagemtimes(u_global, v_global), 1, 1, 1, n);

phy_gs = repmat(gs1_pt / n, n, 1)';
phy_gs = phy_gs + ones(gs1_n, 1) * (0:1 / n:1 - 1 / 2 / n);
phy_gs = repmat(reshape(phy_gs, 1, 1, gs1_n, n), 2, 2, 1, 1);

r_global = bd_r(T, phy_gs);

ia = reshape(repmat(s(1:2, bound_idx), 2, 1), 2, 2 * n);
ja = repmat(reshape(s(1:2, bound_idx), 1, 2 * n), 2, 1);
va = coeff_c * tensorprod(r_global .* uv_global, gs1_wt, 3, 2) / n;
va = reshape(va, 2, 2 * n);

C = C + sparse(ia, ja, va, dof, dof);

v_global = repmat([1 - gs1_pt; gs1_pt], n, 1);

py = p(2, s(1:2, bound_idx));

phy_gs = repmat(gs1_pt / n, 2 * n, 1);
phy_gs = phy_gs + [0; py(1:end - 1)'] * ones(1, gs1_n);

q_global = bd_q(T, phy_gs);

Tr = sparse(py * n + 1, 1:2 * n, ones(1, 2 * n), n + 1, 2 * n);

loc = [s(1, bound_idx) s(2, end)];

F(loc) = F(loc) + coeff_c * Tr * (q_global .* v_global) * gs1_wt' / n;