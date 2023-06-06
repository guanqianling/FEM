%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are six programs in this code package. There are designed to solve Poisson equation and Heat equation with mixed boundary condition.

For different equations, the main program is shown below.

1-D Poisson equation: linear finite element:    FEM_1d_linear_Poisson.m
		      quadratic finite element: FEM_1d_quadratic_Poisson.m

2-D Poisson equation: linear finite element:    FEM_2d_linear_Poisson\FEM_2d_linear_Poisson.m
		      quadratic finite element: FEM_2d_quadratic_Poisson\FEM_2d_quadratic_Poisson.m

2-D Heat equation:    linear finite element:    FEM_2d_linear_Heat\FEM_2d_linear_Heat.m
		      quadratic finite element: FEM_2d_quadratic_Heat\FEM_2d_quadratic_Heat.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Taking 2-D Poisson equation for example, the main program is FEM_2d_linear_Poisson.m. The syntax is 

	z = FEM_2d_linear_Poisson(n);

The input argument(s):

    n is the mesh grid parameter
    (T_end is the final time for Heat equation)

The output argument:

    z is a 1*3 vector, containing the error in L2 norm, L_infty norm and H1 semi-norm

In addition, if one wants obtain the errors and orders of convergence for different mesh sizes, please run order.m after taking an appropriate argument.
If one wants to visualize grid information in 2D problem, please change the variable draw=1 in uniform_mesh.m.
If one wants to visualize numerical solution in 2D problem, please change the varible draw=1 in pre_setting.m.
If one wants to change the Gauss integration, please change it in pre_setting.m.

There is a complex example in 2-D problem, if one wants to test it, please change it in pre_setting.m.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@Author: Haoning Dang
@Time: 05/31/2023


