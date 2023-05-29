function boundary_condition_N

global p s
global F
global coeff_c bd_p
global gs1_pt gs1_wt gs1_n

bound_idx = find(s(4, :) == -2);
n = length(bound_idx);

v_global = repmat([1 - gs1_pt; gs1_pt], n, 1);

py = p(2, s(1:2, bound_idx));

phy_gs = repmat(gs1_pt / n, 2 * n, 1);
phy_gs = phy_gs + [0; py(1:end - 1)'] * ones(1, gs1_n);

p_global = bd_p(phy_gs);

T = sparse(py * n + 1, 1:2 * n, ones(1, 2 * n), n + 1, 2 * n);

loc = [s(1, bound_idx) s(2, end)];

F(loc) = F(loc) + coeff_c * T * (p_global .* v_global) * gs1_wt' / n;

