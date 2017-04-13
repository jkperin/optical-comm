clear, clc, close all

addpath ../f


Fiber = fiber(10e3);

alpha = 2;

% First zero in the frequency response of the fiber occurs when

theta_pos = -2*atan(alpha - sqrt(alpha^2+1));
theta_neg = -2*atan(alpha + sqrt(alpha^2+1));

fband = 25e9;
lamb = ((-200:100) + 1310)*1e-9;
L = [5 10 15 20]*1e3;

figure, hold on, box on
leg = {};
for k = 1:length(L)
    Fiber.L = L(k);

    theta = -1/2*Fiber.beta2(lamb)*(2*pi*fband)^2*Fiber.L;

    plot(lamb*1e9, theta)
    leg = [leg sprintf('L = %d km', Fiber.L/1e3)];
end

legend(leg, 'Location', 'SouthEast')
plot([lamb(1) lamb(end)]*1e9, theta_pos*[1 1], 'k')
plot([lamb(1) lamb(end)]*1e9, theta_neg*[1 1], 'k')
axis([lamb([1 end])*1e9 -3 1])
xlabel('\lambda (nm)')
ylabel('\theta')
