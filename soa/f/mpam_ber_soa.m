function Pe = mpam_ber_soa(mpam, soa, rx, G, PrecdBm)

if nargin == 5
    tx.Prec = 1e-3*10^(PrecdBm/10);
    
    mpam.b = mpam.b*tx.Prec/mean(mpam.a);
    mpam.a = mpam.a*tx.Prec/mean(mpam.a);
    
    soa.G = G;
end    

x = 10*soa.G*mpam.a(end)*linspace(-1, 1, 2^14);
 
fEsp = 1/(soa.Seq(soa.G)*mpam.Rs/2)*pdf('chi2', x/(soa.Seq(soa.G)*mpam.Rs/2), 2);

%% Level 1
% x = 5*sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(1)*soa.Seq(soa.G)*mpam.Rs)*linspace(-1, 1, 2^12);

% fEsp = pdf('chi2', x/(soa.Seq(soa.G)*mpam.Rs/2), 2);

fGauss = pdf('norm', x, soa.G*mpam.a(1), sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(1)*soa.Seq(soa.G)*mpam.Rs));

fx = conv(fEsp, fGauss, 'same');
fx = fx/trapz(x, fx);

% fx = fGauss;
% fisoa.Gure, hold on
% plot(x, fx)

% assert(trapz(x, fx) > 0.9, 'ransoa.Ge for pdf is not enousoa.Gh: %.5f\n', trapz(x, fx))

xe = (x > soa.G*(mpam.b(1)));

Pe = trapz(x(xe), fx(xe));


%% Inner levels

for k = 2:mpam.M-1    
%     x = 5*sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(k)*soa.Seq(soa.G)*mpam.Rs)*linspace(-1, 1, 2^12);

%     fEsp = pdf('chi2', x/(soa.Seq(soa.G)*mpam.Rs/2), 2);
    
    fGauss = pdf('norm', x, soa.G*mpam.a(k), sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(k)*soa.Seq(soa.G)*mpam.Rs));
    
    fx = conv(fEsp, fGauss, 'same');
    fx = fx/trapz(x, fx);
%     fx = fGauss;
%     plot(x, fx)
    
%     assert(trapz(x, fx) > 0.9, 'ransoa.Ge for pdf is not enousoa.Gh: %.5f\n', trapz(x, fx))
    
    xe1 = (x > soa.G*(mpam.b(k)));
    xe2 = (x < soa.G*(mpam.b(k-1)));

    Pe = Pe + trapz(x(xe1), fx(xe1)) + trapz(x(xe2), fx(xe2));   
end

%% last level
% x = 5*sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(end)*soa.Seq(soa.G)*mpam.Rs)*linspace(-1, 1, 2^12);

% fEsp = pdf('chi2', x/(soa.Seq(soa.G)*mpam.Rs/2), 2);
% 
fGauss = pdf('norm', x, soa.G*mpam.a(end), sqrt(rx.N0*mpam.Rs + 2*soa.G*mpam.a(end)*soa.Seq(soa.G)*mpam.Rs));

fx = conv(fEsp, fGauss, 'same');
fx = fx/trapz(x, fx);
% fx = fGauss;

% plot(x, fx)

% assert(trapz(x, fx) > 0.9, 'ransoa.Ge for pdf is not enousoa.Gh: %.5f\n', trapz(x, fx))

xe = (x < soa.G*(mpam.b(end)));

Pe = Pe + trapz(x(xe), fx(xe));

%% Approximate ber
Pe = Pe/mpam.M*1/log2(mpam.M);

