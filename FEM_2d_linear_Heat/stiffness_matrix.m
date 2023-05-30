function stiffness_matrix

global dof_loc dof M_pre
global NT
global p t
global K J
global coeff_c
global A M
global Delta_t

P = reshape(p(:, t), 2, dof_loc, NT);
E = P(:, [3 1 2], :) - P(:, [2 3 1], :);
T = [-E(2, :, :); E(1, :, :)];

A_loc = coeff_c * pagemtimes(permute(T, [2, 1, 3]), T) * K / J ^ 2;

i = repmat(reshape(t, 1, dof_loc * NT), dof_loc, 1);
j = reshape(repmat(t, dof_loc, 1), dof_loc, dof_loc * NT);
va = reshape(A_loc, dof_loc, dof_loc * NT);
vm = repmat(M_pre * J / Delta_t, 1, NT);

A = sparse(i, j, va, dof, dof);
M = sparse(i, j, vm, dof, dof);