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

One can run code directly, default input arguments are provided in each main program.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

For Heat equation, backward Euler scheme is used to iterate time and the time step \Delta t = 4/h^2 for linear element and \Delta t=8/h^3 for quadratic element. If one wants to change the time step, please change it in partition.m.

In addition, if one wants to obtain the errors and orders of convergence for different mesh sizes, please run order.m after taking an appropriate argument.

If one wants to visualize grid information in 2D problem, please change the variable draw=1 in uniform_mesh.m.
If one wants to visualize numerical solution in 2D problem, please change the varible draw=1 in pre_setting.m.
If one wants to change the Gauss integration, please change it in pre_setting.m.

There is a more complex example in 2-D problem, if one wants to test it, please uncomment it in pre_setting.m. What‘s more, one can write a new example, please create a new subroutine and change the corresponding function handle in pre_setting.m, don't have to make changes in each original subroutine.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Run "FEM_2d_linear_Poisson\FEM_2d_linear_Poisson.m", the default input is n = 32 and the output is 

ans =

     8.616744703416944e-04     3.255923574469577e-03     9.112062874794304e-02

they are the error in L2 norm, L_infty norm and H1 semi-norm.


Run "FEM_2d_linear_Poisson\order.m", default N = 10, which means that running "FEM_2d_linear_Poisson.m" 10 times and the inputs are 2^1=2, 2^2=4, 2^3=8, ... , 2^10=1024 respectively. The output is 

ans =

   0.214754357374225                   0   0.522697377282591                   0   1.425299482625889                   0
   0.054713050397443   1.972730506854645   0.168719660758964   1.631347820912727   0.724209571410286   0.976785939672522
   0.013759489571624   1.991458039880079   0.047841159351298   1.818303843052642   0.363879691608069   0.992945708548744
   0.003445312558493   1.997720176987956   0.012685083520588   1.915119332921098   0.182179664614831   0.998099604941914
   0.000861674470342   1.999420024538764   0.003255923574470   1.961994277292620   0.091120628747943   0.999512323582392
   0.000215440356568   1.999854417687535   0.000824885539418   1.980800985028775   0.045564198192371   0.999877021771869
   0.000053861448937   1.999963577013461   0.000207564635566   1.990633284584118   0.022782585925911   0.999969171439204
   0.000013465447288   1.999990887326892   0.000052050076577   1.995588470081017   0.011391353867906   0.999992286473992
   0.000003366367118   1.999997730140692   0.000013032259916   1.997813003636556   0.005695684548956   0.999998071147739
   0.000000841591923   1.999999754629460   0.000003260486136   1.998928291980128   0.002847843226421   0.999999517752991

Column 1, column 3, column 5 are the errors in L2 norm, L_infty norm and H1 semi-norm for different mesh sizes and column 2, column 4, column 6 are the orders of convergence corresponding errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This code is available at https://github.com/guanqianling/FEM

@Author: Haoning Dang
@Date: 05/31/2023
@Email: haoningdang.xjtu@stu.xjtu.edu.cn


