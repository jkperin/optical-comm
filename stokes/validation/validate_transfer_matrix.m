clear, clc, close all

phi = rand(1, 3)*2*pi;
U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];
U = U1*U2*U3; % = [a -b; b^* a^*], for a and b complex

a = U(1, 1);
b = -U(1, 2);


M = [abs(a)^2, abs(b)^2, -2*real(a*b'), 2*imag(a*b');...
    abs(b)^2, abs(a)^2,  2*real(a*b'), -2*imag(a*b');...
    real(a*b), -real(a*b),  real(a^2) - real(b^2), -imag(a^2) - imag(b^2);...
    imag(a*b), -imag(a*b), imag(a^2) - imag(b^2), real(a^2) + real(b^2)];


Etx = 2*randn(2, 1) + 1j*2*randn(2, 1);
Erx = U*Etx;

Ytx = [abs(Etx(1))^2; abs(Etx(2))^2; real(Etx(1)*Etx(2)'); imag(Etx(1)*Etx(2)')]
Yrx = [abs(Erx(1))^2; abs(Erx(2))^2; real(Erx(1)*Erx(2)'); imag(Erx(1)*Erx(2)')]

M*Ytx

norm(M*Ytx-Yrx)
