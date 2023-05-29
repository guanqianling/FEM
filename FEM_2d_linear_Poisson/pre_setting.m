function pre_setting

global draw
global solu rhs coeff_c diff_u
global bd_p bd_r bd_q bd_g
global dof_loc
global gs1_wt gs1_pt gs1_n

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

dof_loc = 3;

gs1_wt = [0.118463442528095, 0.239314335249683,...
          0.284444444444445, 0.239314335249683, 0.118463442528095];
gs1_pt = [0.046910077030668, 0.230765344947158, 0.5, ...
          0.769234655052842, 0.953089922969332];
gs1_n = 5;

integral_transform(2);