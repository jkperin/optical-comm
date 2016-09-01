function eyediagram(x, n, Ntraces)
if not(exist('Ntraces', 'var'))
    Ntraces = 500;
end

N = length(x);
M = min(Ntraces, floor(N/n));
sel = 1:M*n;

x = reshape(x(sel), n, M);
plot(1:n, x)

a = axis;
axis([1 n a(3) a(4)])