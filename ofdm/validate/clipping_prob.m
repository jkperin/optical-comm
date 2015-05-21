clear, clc, close all

r = 0:0.1:5;

x = randn(2^14,1);
for k = 1:length(r)
    Pc1(k) = mean(x < -r(k));
end

x = randn(2^15,1);
for k = 1:length(r)
    Pc2(k) = mean(x < -r(k));
end
   
x = randn(2^16,1);
for k = 1:length(r)
    Pc3(k) = mean(x < -r(k));
end

plot(r, log10(qfunc(r)), 'k', r, log10(Pc1), 'r', r, log10(Pc2), 'b', r, log10(Pc3), 'c')
xlabel('Clipping ratio (r)')
ylabel('Clipping probability')
legend('Theory', '2^{14}', '2^{15}', '2^{16}')