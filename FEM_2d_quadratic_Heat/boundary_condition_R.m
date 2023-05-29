function boundary_condition_R(T)

global p s
global C F
global dof
global gs1_pre
global coeff_c bd_r bd_q
global gs1_pt gs1_wt gs1_n

bound_idx = find(s(4, :) == -3);
n = length(bound_idx);

u_global = reshape(gs1_pre, 3, 1, gs1_n);
v_global = reshape(gs1_pre, 1, 3, gs1_n);

uv_global = repmat(pagemtimes(u_global, v_global), 1, 1, 1, n);

phy_gs = repmat(gs1_pt / n, n, 1)';
phy_gs = phy_gs + ones(gs1_n, 1) * (0:1 / n:1 - 1 / 2 / n);
phy_gs = repmat(reshape(phy_gs, 1, 1, gs1_n, n), 3, 3, 1, 1);

r_global = bd_r(T, phy_gs);

loc = [s(1, bound_idx); bound_idx + size(p,2)-size(s,2); s(2, bound_idx)];

ia = reshape(repmat(loc, 3, 1), 3, 3 * n);
ja = repmat(reshape(loc, 1, 3 * n), 3, 1);
va = coeff_c * tensorprod(r_global .* uv_global, gs1_wt, 3, 2) / n;
va = reshape(va, 3, 3 * n);

C = C + sparse(ia, ja, va, dof, dof);

v_global = repmat(gs1_pre, n, 1);

Tr = reshape(repmat(0:1 / n:1 - 1 / n, 3, 1), 3 * n, 1);
phy_gs = repmat(gs1_pt / n, 3 * n, 1) + Tr * ones(1, gs1_n);

q_global = bd_q(T, phy_gs);

Tr = sparse(sort([1:2*n + 1 3:2:2*n]), 1:3*n, ones(1, 3*n), 2*n + 1, 3*n);

loc = zeros(1 , 2 * n + 1);
loc(1:2:end) = [s(1, bound_idx) s(2, bound_idx(end))];
loc(2:2:end) = bound_idx + size(p, 2) - size(s, 2);

F(loc) = F(loc) + coeff_c * Tr * (q_global .* v_global) * gs1_wt' / n;