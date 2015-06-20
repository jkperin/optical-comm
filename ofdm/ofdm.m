classdef ofdm 
    properties
        N  % number of subcarriers
        CS % nominal constellation size
        Rs
    end
    
    properties(Constant, GetAccess=private)
        max_iterations = 100; % maximum number of iterations
        frac_incl = 0.9999;    % fraction of energy to be included within CP length
    end
    
    methods
        function ofdm()
        end
        function power_allocation()
            %% Preemphasis at the transmitter
            switch sim.type
                case 'preemphasis'
                    [Pn, CS, aux] = preemphasis(this, tx, fiber, rx, Gch, sim);
                %% Optimal power allocation and bit loading   
                case 'palloc'
                    [Pn, CS, aux] = palloc(this, tx, fiber, rx, Gch, sim);
                otherwise 
                    error('power_allocation: invalid option') 
            end
        end
        
        
        %% Calculate the cyclic prefix length after oversampling 
        % The channel is assumed to be a 2nd-order filter with damping ratio = 1.

        % Algorithm:
        % 1. Assumes a cyclic prefix length k
        % 2. From k, it calculates the new sampling rate fs
        % 3. Calculate the number of samples at both sides (Nneg and Npos) that 
        % contains the desired fraction of energy
        % 4. From that calculates Npre_os = Npos + Neg; (zero is not included)
        % 5. If Npre_os == k end simulation, otherwise increment k and repeat

        function [Npre_os, Nneg, Npos] = cyclic_prefix(this, tx, rx, sim)
        
            Mct = sim.Mct;
            Ntot = 1024;

            Nd = ceil(this.N/this.ros);
            
            gdac = impz(tx.filter.den, tx.filter.num) 
            gadc = impz(rx.filter.den, rx.filter.num)

            k = 0;
            Npre_os = Inf;
            while Npre_os > k && k < this.max_iterations
                fs = this.Rs*(this.Nc + k)/Nd;
                fsct = Mct*fs;
                dt = 1/fsct;

                % Group delay of modulator in samples
                hl_delay = tx.hl_delay*fsct;

                % Channel impulse response
                tct = (0:Ntot-1)*dt;
                hct = tx.hl(tct);
                
                % Total impulse response in continuous time
                pct = conv(conv(gdac, hct, 'full'), gadc, 'full');
                pct = pct/max(pct);

                % Remove group delay due to filters so that impulse response is
                % centered at zero.
                tct = (0:length(pct)-1)*dt;
                tct = tct - ceil(tx.gdac_delay + rx.gadc_delay + hl_delay)*dt;

                n0 = ceil(tx.gdac_delay + hl_delay + rx.gadc_delay) + 1; % new zero after remove group delay

                % Sampling at the chip rate
                t = tct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);
                p = pct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);

                % CP based on energy
                en_frac = cumsum(p.^2)/sum(p.^2);
                Nneg = sum(en_frac >= (1 - frac_incl)/2 & t < 0);
                Npos = sum(en_frac <= (1 + frac_incl)/2 & t > 0);

                % Number of samples necessary to attain desired fraction of energy
                Npre_os = Nneg + Npos;

                k = k + 1;
            end

            assert(k ~= max_it, 'CP calculation did not converge');
        end
    end
end