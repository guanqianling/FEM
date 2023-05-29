function draw_uh

global NT
global p t
global uh

figure
hold on
for i = 1 : NT
    trimesh([1 2 3], p(1, t(1:3, i)), p(2, t(:, i)), uh(t(:, i)));
end
hold off
view(-38,30)
