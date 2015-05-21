% Gain distribution
% n = primary, r = secondary
% pnnr = p(n|n + r)
function pnnr = apd_gain_distribution(n, keff, Gapd)
% 
% prv = @(v, r) v*((1-keff)*(Gapd-1)/Gapd).^(r-v).*gamma(r/(1-keff)).*((1 + keff*(Gapd-1))/Gapd).^((v+keff*(r-v))/(1-keff))./...
%     ((v + keff*(r-v)).*factorial(r-v).*((v + keff*(r-v))/(1-keff)));


T1 = (n-1)*log((1-keff)*(Gapd-1)/Gapd);

T2 = log(gamma(n/(1-keff)));

T3 = ((1 + keff*(n-1))/(1-keff)).*log((1 + keff*(Gapd-1))/Gapd);

T4 = log((1 + keff*(n-1)).*factorial(n-1));
% T4 = log((1 + keff*(n-1))) + sum(log(1:n-1));

T5 = log(gamma((1 + keff*(n-1))/(1-keff)));

pnnr = T1 + T2 + T3 - T4 - T5;

pnnr = exp(pnnr);

%% Generic equation with Stirling's approximation for the factorial and gamma functions
% Approximating the factorial and Gamma functions
% 
% r = n;
% n = 1;
% f = (n+r)./(n*Gapd);
% X = f - 1;
% 
% T1 = log(n./sqrt(2*pi*(n + r).*(n + keff*r).*r));
% T2 = r.*log((1 - n.*X./r));
% T3 = ((n + keff*r)./(1 - keff)).*log(1 + n.*X*(1-keff)./(n + keff*r));
% 
% pnnr = T1 + T2 + T3;
% 
% % if any(isinf(pnnr) | isnan(pnnr)) 
% %     error('gainpmf', 'Error while calculating gain pmf')
% % end
% 
% pnnr = exp(pnnr);
% 
