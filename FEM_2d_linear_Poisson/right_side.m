function right_side

global F
global dof dof_loc
global p t NT
global gs2_pt gs2_wt
global phy_x phy_y
global u_global
global rhs
global K

px = p(1, :);
py = p(2, :);

area_xyz = [1 - sum(gs2_pt); gs2_pt];
phy_x = px(t)' * area_xyz;
phy_y = py(t)' * area_xyz;

idx = reshape(ones(dof_loc, 1) * (1:NT), dof_loc * NT, 1);

F_global = rhs(phy_x, phy_y);
F_global = F_global(idx, :);

gs2_ref = [1 - gs2_pt(1,:)' - gs2_pt(2,:)', gs2_pt(1,:)', gs2_pt(2,:)']';
u_global = repmat(gs2_ref, NT, 1);

loc = reshape(t, 1, 3 * size(t, 2));
T = sparse(loc, 1:dof_loc * NT, ones(1, dof_loc * NT), dof, dof_loc * NT);

F = T * (F_global .* u_global) * gs2_wt' .* K;