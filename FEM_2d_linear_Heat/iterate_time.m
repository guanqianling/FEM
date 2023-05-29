function iterate_time

global uh
global u0
global p
global T Delta_t

uh = u0(p(1, :)', p(2, :)');

N = T / Delta_t;

for i = 1 : N
    ti = i * Delta_t;
    
    right_side(ti);

    boundary_condition_R(ti);
    
    boundary_condition_N(ti); 

    boundary_condition_D(ti);

    solve_CF;
end




