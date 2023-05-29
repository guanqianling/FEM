function draw_uh

global NT
global p t
global uh

figure
hold on
for i = 1 : NT
    trimesh([1 2 3], p(1, t(1:3, i)), p(2, t(1:3, i)), uh(t(1:3, i)));
%     patch('faces',[1 6 2 4 3 5],'vertices',[p(1, t(:, i))' p(2, t(:, i))' uh(t(:, i))],'FaceColor','none');
end
hold off
view(-38,30)
