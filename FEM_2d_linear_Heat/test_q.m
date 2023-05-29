function z = test_q(t, y)
z = exp(t + 1 + 2 * y) + t * cos(1) + t .^ 2 .* y .^ 2 .* (exp(t + 1 + 2 * y) + t * sin(1) + t * cos(2 * y));