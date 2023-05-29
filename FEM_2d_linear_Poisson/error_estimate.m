function error_estimate

global draw
global solu diff_u
global phy_x phy_y
global dof_loc
global NT
global p t
global u_global
global uh
global gs2_wt gs2_n
global K
global err_L2 err_Linf err_H1_semi

if draw == 1
    draw_uh;
end

u_real_global = solu(phy_x, phy_y);

loc = reshape(t, 3 * size(t, 2), 1);

row = reshape(repmat(1:NT, dof_loc, 1), 1, dof_loc * NT);
Tr = sparse(row, 1:dof_loc * NT, ones(1, dof_loc * NT), NT, dof_loc * NT);

u_global = Tr * (repmat(uh(loc), 1, gs2_n) .* u_global);

err = u_real_global - u_global;
err_L2 = sqrt(sum(err .* err * gs2_wt' .* K));

err_Linf = max(max(abs(err)));

diff_u_real_global = diff_u(phy_x, phy_y);

ia = repmat(1:dof_loc * NT, dof_loc, 1);
ja = reshape(repmat(1:dof_loc:dof_loc * NT, dof_loc, 1), 1, dof_loc * NT);
ja = repmat(ja, dof_loc, 1) + repmat(0:dof_loc - 1, dof_loc * NT, 1)';
va = [ones(1, dof_loc * NT); p(:, reshape(t, dof_loc * NT, 1))];

M = sparse(ia, ja, va, dof_loc * NT, dof_loc * NT);
b = uh(reshape(t, dof_loc * NT, 1));
x = M \ b;

diff_u_global_x = repmat(x(2:dof_loc:dof_loc * NT), 1, gs2_n);
diff_u_global_y = repmat(x(3:dof_loc:dof_loc * NT), 1, gs2_n);

diff_u_global = [diff_u_global_x; diff_u_global_y];

err = diff_u_real_global - diff_u_global;
err_H1_semi = sqrt(sum(err .* err * gs2_wt' * K));