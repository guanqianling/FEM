function partition(n)

global p s t
global Delta_t
global NT
global K
global dof

[p, s, t] = uniform_mesh(n);

Delta_t = 8 / n ^ 3;

NT = 2 * n * n;

K = 1 / 2 / n / n;

dof = (2 * n + 1) ^ 2;