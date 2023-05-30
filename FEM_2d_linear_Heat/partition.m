function partition(n)

global p s t
global Delta_t
global NT
global K J
global dof

[p, s, t] = uniform_mesh(n);

Delta_t = 4 / n ^ 2;

NT = 2 * n * n;

K = 1 / 2 / n / n;

J = 1 / n / n;

dof = (n + 1) * (n + 1);