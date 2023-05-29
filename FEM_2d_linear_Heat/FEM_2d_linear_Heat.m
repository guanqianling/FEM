function z = FEM_2d_linear_Heat(n, T_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is designed to solve Heat equation with mixed boundary
% condition in $[0,T]\times\Omega$ and $\Omega=[0,1]\times[0,1]$
% u_t-\nabla \cdot(c\nabla u)=f in [0,T]\times\Omega,
%                    u(0,x,y)=u_0 at t=0 and in \Omega,
%    \nabla u\cdot \mathbf{n}=p on [0,T]\times \{x=0\},
% \nabla u\cdot \mathbf{n}+ru=q on [0,T]\times \{x=1\},
%                           u=g on [0,T]\times \{y=0,1\}
% 
%----The input arguments:
% n is the mesh grid parameter, default = 32
% T_end is the final time, default = 1
%----The output arguments:
% z(1) is the error corresponding to the L2 norm
% z(2) is the error corresponding to the L_infty norm
% z(3) is the error corresponding to the H1 semi-norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    T_end = 1;
end
if nargin < 1
    n = 32;
end

clear global

global err_L2 err_Linf err_H1_semi
global T

T = T_end;

pre_setting;

tic

partition(n);

stiffness_matrix;

iterate_time;

toc

error_estimate;

z = [err_L2 err_Linf err_H1_semi];

clear global