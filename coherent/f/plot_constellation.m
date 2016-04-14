function plot_constellation(Y, dataTX, M)

box on, hold on
for k = 0:M-1
    s = Y(dataTX == k);
    scatter(real(s), imag(s), 10, 'o', 'filled')
end

ref = qammod(0:M-1, M);
scatter(real(ref), imag(ref), 20, 'x')

axis([-2 2 -2 2]);
axis square