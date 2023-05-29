function z = test_f(t, x, y)
z = -4 * exp(t + x + 2 * y) + (t + 1) .* sin(x) + (4 * t + 1) .* cos(2 * y);