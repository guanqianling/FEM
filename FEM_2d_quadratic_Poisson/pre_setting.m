function pre_setting

global draw
global solu rhs coeff_c diff_u
global bd_p bd_r bd_q bd_g
global dof_loc A_pre gs1_pre gs2_pre
global gs1_wt gs1_pt gs1_n
global gs2_pt

format long e

draw = 0;

coeff_c = 1;

solu = @u_real;
rhs = @f;
diff_u = @grad_u_real;
bd_p = @p;
bd_r = @r;
bd_q = @q;
bd_g = @g;

% solu = @test_u_real;
% rhs = @test_f;
% diff_u = @test_grad_u_real;
% bd_p = @test_p;
% bd_r = @test_r;
% bd_q = @test_q;
% bd_g = @test_g;

dof_loc = 6;
A_pre(:, :, 1) = [ 1/2,  1/6, 0,    0,    0, -2/3;
                   1/6,  1/2, 0,    0,    0, -2/3;
                     0,    0, 0,    0,    0,    0;
                     0,    0, 0,  4/3, -4/3,    0;
                     0,    0, 0, -4/3,  4/3,    0;
                  -2/3, -2/3, 0,    0,    0,  4/3];
A_pre(:, :, 2) = [ 1/2, 0,  1/6,    0, -2/3,    0;
                   1/6, 0, -1/6,  2/3,    0, -2/3;
                     0, 0,    0,    0,    0,    0;
                     0, 0,  2/3,  2/3, -2/3, -2/3;
                     0, 0, -2/3, -2/3,  2/3,  2/3;
                  -2/3, 0,    0, -2/3,  2/3,  2/3];
A_pre(:, :, 3) = [ 1/2,  1/6, 0,    0,    0, -2/3;
                     0,    0, 0,    0,    0,    0;
                   1/6, -1/6, 0,  2/3, -2/3,    0;
                     0,  2/3, 0,  2/3, -2/3, -2/3;
                  -2/3,    0, 0, -2/3,  2/3,  2/3;
                     0, -2/3, 0, -2/3,  2/3,  2/3];
A_pre(:, :, 4) = [ 1/2, 0,  1/6,    0, -2/3,    0;
                     0, 0,    0,    0,    0,    0;
                   1/6, 0,  1/2,    0, -2/3,    0;
                     0, 0,    0,  4/3,    0, -4/3;
                  -2/3, 0, -2/3,    0,  4/3,    0;
                     0, 0,    0, -4/3,    0,  4/3];

gs1_wt = [0.118463442528095, 0.239314335249683,...
          0.284444444444445, 0.239314335249683, 0.118463442528095];
gs1_pt = [0.046910077030668, 0.230765344947158, 0.5, ...
          0.769234655052842, 0.953089922969332];
gs1_n = 5;
integral_transform(2);

gs1_pre = [2 * gs1_pt .^ 2 - 3 * gs1_pt + 1; 
           4 * gs1_pt .* (1 - gs1_pt); 
           gs1_pt .* (2 * gs1_pt - 1)];
gs2_pre = [2 * gs2_pt(1,:)' .* gs2_pt(1,:)' + ...
               2 * gs2_pt(2,:)' .* gs2_pt(2,:)' + ...
               4 * gs2_pt(1,:)' .* gs2_pt(2,:)' - ...
               3 * gs2_pt(1,:)' - ...
               3 * gs2_pt(2,:)' + ...
               1, ...
           gs2_pt(1,:)' .* (2 * gs2_pt(1,:)' - 1), ...
           gs2_pt(2,:)' .* (2 * gs2_pt(2,:)' - 1), ...
           4 * gs2_pt(1,:)' .* gs2_pt(2,:)', ...
           4 * gs2_pt(2,:)' .* (1 - gs2_pt(1,:)' - gs2_pt(2,:)'), ...
           4 * gs2_pt(1,:)' .* (1 - gs2_pt(1,:)' - gs2_pt(2,:)')]';