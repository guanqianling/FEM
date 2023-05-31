function res = order

N = 8;
n = 2 .^ (1:N);

err = zeros(N, 3);

for i = 1 : N
%     a = FEM_1d_linear_Poisson(n(i));
    a = FEM_1d_quadratic_Poisson(n(i));
    err(i, :) = a(1:3);
end

format long

err_reduce = abs(err(2:N, :) ./ err(1:N-1, :));

h_reduce = (n(1:N - 1) ./ n(2:N))' * ones(1, 3);

err_order = log(err_reduce) ./ log(h_reduce);

res = [err(:, 1) [0; err_order(:, 1)] err(:, 2) [0; err_order(:, 2)] ...
       err(:, 3) [0; err_order(:, 3)]];

errorder.data = res;
errorder.tableColLabels = {'$\|e_h\|_0$','order','$\|e_h\|_\infty$','order','$|e_h|_1$','order'};
errorder.tableRowLabels = {'1/2','1/4','1/8','1/16','1/32','1/64','1/128','1/256'};
errorder.dataFormat = {'%.4e',1,'%.4f',1,'%.4e',1,'%.4f',1,'%.4e',1,'%.4f',1};
latexTable(errorder);