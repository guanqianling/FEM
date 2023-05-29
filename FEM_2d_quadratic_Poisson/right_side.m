function right_side

global F
global dof dof_loc
global p t NT
global gs2_pt gs2_wt
global phy_x phy_y
global gs2_pre
global u_global
global rhs
global K

px = p(1, :);
py = p(2, :);

area_xyz = [1 - sum(gs2_pt); gs2_pt];
phy_x = px(t(1:3, :))' * area_xyz;
phy_y = py(t(1:3, :))' * area_xyz;

idx = reshape(ones(dof_loc, 1) * (1:NT), dof_loc * NT, 1);

F_global = rhs(phy_x, phy_y);
F_global = F_global(idx, :);

u_global = repmat(gs2_pre, NT, 1);

loc = reshape(t, 1, dof_loc * size(t, 2));
T = sparse(loc, 1:dof_loc * NT, ones(1, dof_loc * NT), dof, dof_loc * NT);

F = T * (F_global .* u_global) * gs2_wt' .* K;