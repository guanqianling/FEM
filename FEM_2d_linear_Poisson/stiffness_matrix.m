function stiffness_matrix

global dof_loc dof
global NT
global p t
global K
global coeff_c
global A

ia = zeros(dof_loc, dof_loc * NT)';
ja = zeros(dof_loc, dof_loc * NT)';
va = zeros(dof_loc, dof_loc * NT)';

for i = 1 : NT
    P = p(:, t(:, i));
    E = P(:, [3 1 2]) - P(:, [2 3 1]);
    
    T = [-E(2, :); E(1, :)];
    A_loc = coeff_c * (T.'*T) * K / det([E(:, 3) -E(:, 2)]) ^ 2;

    [x, y] = meshgrid(t(:, i), t(:, i));
    ia(dof_loc * (i - 1) + (1:dof_loc), :) = x;
    ja(dof_loc * (i - 1) + (1:dof_loc), :) = y;
    va(dof_loc * (i - 1) + (1:dof_loc), :) = A_loc;
end

A = sparse(ia, ja, va, dof, dof);