function z = test_q(y)
z = exp(1 + 2 * y) + cos(1) + y .^ 2 .* (exp(1 + 2 * y) + sin(1) + cos(2 * y));