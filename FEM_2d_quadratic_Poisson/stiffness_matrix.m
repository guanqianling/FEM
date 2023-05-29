function stiffness_matrix

global dof_loc dof A_pre
global NT
global p t
global coeff_c K
global A

ia = zeros(dof_loc, dof_loc * NT)';
ja = zeros(dof_loc, dof_loc * NT)';
va = zeros(dof_loc, dof_loc * NT)';

for i = 1 : NT
    P = p(:, t(1:3, i));
    E = P(:, [3 1 2]) - P(:, [2 3 1]);

    D = inv([E(:, 3) -E(:, 2)])';
    T = D' * D * K * 2;

    A_loc = T(1, 1) * A_pre(:, :, 1) + T(1, 2) * A_pre(:, :, 2) + ...
            T(2, 1) * A_pre(:, :, 3) + T(2, 2) * A_pre(:, :, 4);

    [x, y] = meshgrid(t(:, i), t(:, i));
    ia(dof_loc * (i - 1) + (1:dof_loc), :) = x;
    ja(dof_loc * (i - 1) + (1:dof_loc), :) = y;
    va(dof_loc * (i - 1) + (1:dof_loc), :) = coeff_c * A_loc;
end

A = sparse(ia, ja, va, dof, dof);