function boundary_condition_D(T)

global p s
global C F
global bd_g

bound_idx = find(s(4, :) == -1);
bound_idx1 = reshape(s(1:2, bound_idx), 1, 2 * length(bound_idx));
bound_idx = [unique(sort(bound_idx1)) bound_idx + size(p, 2) - size(s, 2)];

C(bound_idx, :) = 0;
C(bound_idx, bound_idx) = eye(length(bound_idx));

F(bound_idx, 1) = bd_g(T, p(1, bound_idx), p(2, bound_idx));