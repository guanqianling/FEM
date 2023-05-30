function stiffness_matrix

global dof_loc dof
global NT
global p t
global K J
global coeff_c
global A

P = reshape(p(:, t), 2, dof_loc, NT);
E = P(:, [3 1 2], :) - P(:, [2 3 1], :);
T = [-E(2, :, :); E(1, :, :)];

A_loc = coeff_c * pagemtimes(permute(T, [2, 1, 3]), T) * K / J ^ 2;

i = repmat(reshape(t, 1, dof_loc * NT), dof_loc, 1);
j = reshape(repmat(t, dof_loc, 1), dof_loc, dof_loc * NT);
v = reshape(A_loc, dof_loc, dof_loc * NT);

A = sparse(i, j, v, dof, dof);