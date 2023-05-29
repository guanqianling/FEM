function boundary_condition_D

global p s
global A F
global bd_g

bound_idx = find(s(4, :) == -1);
bound_idx = reshape(s(1:2, bound_idx), 1, 2 * length(bound_idx));
bound_idx = unique(sort(bound_idx));

A(bound_idx, :) = 0;
A(bound_idx, bound_idx) = eye(length(bound_idx));

F(bound_idx, 1) = bd_g(p(1, bound_idx), p(2, bound_idx));