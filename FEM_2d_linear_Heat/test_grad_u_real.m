function z = test_grad_u_real(t, x, y)
z = [exp(t + x + 2 * y) + t .* cos(x); 2 * exp(t + x + 2 * y) - 2 * t .* sin(2 * y)];