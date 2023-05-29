function z = test_g(t, x, y)
z = exp(t + x + 2 * y) + t .* sin(x) + t .* cos(2 * y);