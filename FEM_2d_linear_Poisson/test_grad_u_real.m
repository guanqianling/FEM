function z = test_grad_u_real(x, y)
z = [exp(x + 2 * y) + cos(x); 2 * exp(x + 2 * y) - 2 * sin(2 * y)];