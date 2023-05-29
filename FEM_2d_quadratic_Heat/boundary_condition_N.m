function boundary_condition_N(T)

global p s
global F
global gs1_pre
global coeff_c bd_p
global gs1_pt gs1_wt gs1_n

bound_idx = find(s(4, :) == -2);
n = length(bound_idx);

v_global = repmat(gs1_pre, n, 1);

Tr = reshape(repmat(0:1 / n:1 - 1 / n, 3, 1), 3 * n, 1);
phy_gs = repmat(gs1_pt / n, 3 * n, 1) + Tr * ones(1, gs1_n);

p_global = bd_p(T, phy_gs);

Tr = sparse(sort([1:2*n + 1 3:2:2*n]), 1:3*n, ones(1, 3*n), 2*n + 1, 3*n);

loc = zeros(1 , 2 * n + 1);
loc(1:2:end) = [s(1, bound_idx) s(2, bound_idx(end))];
loc(2:2:end) = bound_idx + size(p, 2) - size(s, 2);

F(loc) = F(loc) + coeff_c * Tr * (p_global .* v_global) * gs1_wt' / n;