%% Calculate OSNR required for a give fiber length
clear, clc, close all

addpath f/ % Juniper project specific functions
addpath ../mpam % PAM
addpath ../f % general functions
addpath ../soa % for pre-amplifier 
addpath ../apd % for PIN photodetectors

BERtarget = 1.8e-4;

Lkm = 0:20;
alpha = -1;
Ntaps = 5;

for k = 1:length(Lkm)
    fprintf('===================== L = %d km =====================\n', Lkm(k))
    [BER, OSNRdB] = pam_sim_fun(Lkm(k), alpha, Ntaps);
    
    data.Lkm(k) = Lkm(k);
    data.BER{k} = BER;
    data.OSNRdB{k} = OSNRdB;
    
    try
        idx1 = find(BER.gauss ~= 0);
        idx2 = find(BER.count ~= 0);

        [x1, ~, exitflag1] = fzero(@(x) interp1(OSNRdB(idx1), log10(BER.gauss(idx1)), x, 'linear', 'extrap') - log10(BERtarget), OSNRdB(idx1(1)));
        if exitflag1 ~= 1
            disp('Gauss approx exit flag different from 1')
            exitflag1
        end

        OSNRdBReq.gauss(k) = x1;

        [x2, ~, exitflag2] = fzero(@(x) interp1(OSNRdB(idx2), log10(BER.count(idx2)), x, 'linear', 'extrap') - log10(BERtarget), OSNRdB(idx2(1)));
        if exitflag2 ~= 1
            disp('Count exit flag different from 1')
            exitflag2
        end

        OSNRdBReq.gauss(k) = x1; 
        OSNRdBReq.count(k) = x2;
    catch e
        disp('An error occurred')
        disp(e.message)
    end
end

figure, hold on, box on
plot(Lkm, OSNRdBReq.gauss, 'LineWidth', 4)
plot(Lkm, OSNRdBReq.count, 'LineWidth', 4)
xlabel('Fiber length (km)', 'FontSize', 18);
ylabel('OSNR required (dB)', 'FontSize', 18);
set(gca, 'FontSize', 18)
legend('Gauss', 'Count')