function OSNRdB = estimate_osnr(E, lambda, f, verbose)

[YdBm, ly] = osa(E, lambda, f, 0.1e-9, verbose);

[SandN, idx] = max(YdBm); % signal + noise power

tolVec = [1e-3 1e-2 1e-1 1];
% try different tolerances
for k = 1:length(tolVec)
    tol = tolVec(k);
    idx1 = find(abs(diff(YdBm(idx:-1:1))) < tol); % Select points where difference was below tol
    idx2 = find(abs(diff(YdBm(idx:1:end))) < tol); % Select points where difference was below tol
    if not(isempty(idx1)) && not(isempty(idx2))
        break
    end
end
        
if isempty(idx1) || isempty(idx2)
    warning('estimate_osnr: it was not possible to find noise floor')
    OSRNdB = [];
    return;
end

idx1 = idx - idx1(1);
idx2 = idx + idx2(1);

NdB = interp1([ly(idx1), ly(idx2)], [YdBm(idx1), YdBm(idx2)], ly(idx));

SN = 10^(SandN/10);
N = 10^(NdB/10);

OSNRdB = 10*log10(SN/N - 1);

if exist('verbose', 'var') && verbose
    figure(132), hold on
    h1 = plot(1e9*ly(idx), SandN, 'ok');
    plot(1e9*[ly(idx1(1)) ly(idx2(1))], [YdBm(idx1(1)), YdBm(idx2(1))], 'x-k')
    plot(1e9*ly(idx), NdB, 'ok')
    legend(h1, sprintf('OSNR = %.2f dB', OSNRdB));
end
