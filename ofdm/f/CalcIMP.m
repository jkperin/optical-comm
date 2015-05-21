%% Calculate intermodulation products

function [IMP, S1] = CalcIMP(Pn, ofdm, tx, fiber, rx, sim)

J = @(n, u) (n == 0).*(1 - u.^2/4) + (abs(n) == 1).*n.*u/2;
% J = @(n, u) (n == 0).*besselj(0, u) + (abs(n) == 1).*((-1).^n).*besselj(1, u);

W = 2*pi*ofdm.fc;
beta2 = fiber.beta2;

% sigtx = sqrt(2*sum(Pn));
sigtx = sqrt(2*sum(Pn.*abs(tx.Gdac).^2));

I0 = rx.R*tx.kappa*(sigtx*sim.rcliptx);
mIM = 2*sqrt(Pn.*abs(tx.Gdac.*tx.kappa.*tx.Hmod).^2)/I0;
mFM = tx.alpha*mIM/2;

A2 = [-1 1; 1 1];
A3 = [-1, -1, +1; -1, +1, +1; +1, -1, +1; +1, +1, -1; +1, +1, +1];
A4 = [-1, -1, -1, +1; -1, -1, +1, +1; -1, +1, -1, +1; -1, +1, +1, +1; 1, -1, -1, +1;
     +1, -1, +1, +1; +1, +1, -1, +1; +1, +1, +1, -1; +1, +1, +1, +1];

S1 = zeros(1, ofdm.Nu/2);
S2 = zeros(1, ofdm.Nu/2);
S3 = zeros(1, ofdm.Nu/2);
S4 = zeros(1, ofdm.Nu/2);
for n1 = 1:ofdm.Nu/2
    n = zeros(1, ofdm.Nu/2);
    n(n1) = 1;
    
    Wimp = sum(W.*n);
    
    theta = -0.5*beta2*W*Wimp*fiber.L;
    
    u = 2*mFM.*sin(theta);
    
    S1(n1) = abs(I0*prod(J(n, u))*(1 - sum(mIM.*cos(theta).*n./u)))^2;
    
    for n2 = n1+1:ofdm.Nu/2 
        for kk = 1:size(A2, 1)
            n = zeros(1, ofdm.Nu/2);
            n([n1 n2]) = A2(kk, :);

            Wimp = sum(W.*n);

            theta = -0.5*beta2*W*Wimp*fiber.L;

            u = 2*mFM.*sin(theta);

            nn = sum((1:ofdm.Nu/2).*n);
            
            if nn <= 0 || nn > ofdm.Nu/2
                continue;
            else
                S2(nn) = S2(nn) + abs(I0*prod(J(n, u))*(1 - sum(mIM.*cos(theta).*n./u))).^2;
            end
        end
        
%         for n3 = n2+1:ofdm.Nu/2
%             for kk = 1:size(A3, 1)
%                 n = zeros(1, ofdm.Nu/2);
%                 n([n1 n2 n3]) = A3(kk, :);
% 
%                 Wimp = sum(W.*n);
% 
%                 theta = -0.5*beta2*W*Wimp*fiber.L;
% 
%                 u = 2*mFM.*sin(theta);
% 
%                 nn = sum((1:ofdm.Nu/2).*n);
%                 
%                 if nn <= 0 || nn > ofdm.Nu/2
%                     continue;
%                 else
%                     S3(nn) = S3(nn) + abs(I0*prod(J(n, u))*(1 - sum(mIM.*cos(theta).*n./u))).^2;
%                 end
%             end
            
%             for n4 = n3+1:ofdm.Nu/2
%                 for kk = 1:size(A4, 1)
%                     n = zeros(1, ofdm.Nu/2);
%                     n([n1 n2 n3 n4]) = A4(kk, :);
% 
%                     Wimp = sum(W.*n);
% 
%                     theta = -0.5*beta2*W*Wimp*fiber.L;
% 
%                     u = 2*mFM.*sin(theta);
% 
%                     nn = sum((1:ofdm.Nu/2).*n);
% 
%                     if nn <= 0 || nn > ofdm.Nu/2
%                         continue;
%                     else
%                         S4(nn) = S4(nn) + abs(I0*prod(J(n, u))*(1 - sum(mIM.*cos(theta).*n./u))).^2;
%                     end
%                 end 
%             end
%         end
    end
end

IMP = (S2 + S3 + S4);
% 
% figure
% stem(1:ofdm.Nu/2, IMP.*abs(rx.Gadc).^2, '*r')



