function stiffness_matrix

global dof_loc dof A_pre
global NT
global p t
global coeff_c K
global A

P = reshape(p(:, t(1:3, :)), 2, 3, NT);
E = P(:, [3 1 2], :) - P(:, [2 3 1], :);

D = pageinv([E(:, 3, :) -E(:, 2, :)]);
T = pagemtimes(D, permute(D, [2, 1, 3])) * K * 2;
T = reshape(T, 4, NT);

i = repmat(reshape(t, 1, dof_loc * NT), dof_loc, 1);
j = reshape(repmat(t, dof_loc, 1), dof_loc, dof_loc * NT);
v = coeff_c * reshape(tensorprod(A_pre, T, 3, 1), dof_loc, dof_loc * NT);

A = sparse(i, j, v, dof, dof);