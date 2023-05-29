function z = FEM_2d_quadratic_Poisson(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is designed to solve Possion equation with mixed boundary
% condition in $\Omega=[0,1]\times[0,1]$
%    -\nabla \cdot(c\nabla u)=f in \Omega,
%    \nabla u\cdot \mathbf{n}=p on x=0,
% \nabla u\cdot \mathbf{n}+ru=q on x=1,
%                           u=g on y=0,1
% 
%----The input arguments:
% n is the mesh grid parameter, default = 32
%----The output arguments:
% z(1) is the error corresponding to the L2 norm
% z(2) is the error corresponding to the L_infty norm
% z(3) is the error corresponding to the H1 semi-norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    n = 32;
end

clear global

global err_L2 err_Linf err_H1_semi

pre_setting;

tic

partition(n);

stiffness_matrix;

right_side;

boundary_condition_R;

boundary_condition_N;

boundary_condition_D;

solve_AF;

toc

error_estimate;

z = [err_L2 err_Linf err_H1_semi];

clear global