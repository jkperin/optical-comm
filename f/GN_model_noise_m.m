function NL = GN_model_noise_m(P, D)
%% Compute nonlinear noise power at each channel
% Inputs:
% - P: input power at each span
% - D: cell array of nonlinear coefficients. D is indexed such that D{l+2},
% where l = -1, 0, 1.
% Output:
% - NL: nonlinear noise power

N = length(P);
NL = zeros(size(P));
Ncenter = (size(D{1}, 1)+1)/2;
parfor n = 1:N
    if P(n) == 0
        continue
    end
    for n1 = 1:N
        if P(n1) == 0
            continue
        end
        for n2 = 1:N
            for l = -1:1
                idx = n1 + n2 - n + l;
                if idx < 1 || idx > N
                    continue
                end                
                if P(n2) == 0 || P(idx) == 0
                    continue
                end
                
                i = Ncenter - (n1 - n);
                j = Ncenter + (n2 - n);
                NL(n) = NL(n) + P(n1)*P(n2)*P(idx)*D{l+2}(i, j);
            end
        end
    end
end
