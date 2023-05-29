function partition(n)

global p s t
global NT
global K
global dof

[p, s, t] = uniform_mesh(n);

NT = 2 * n * n;

K = 1 / 2 / n / n;

dof = (n + 1) * (n + 1);