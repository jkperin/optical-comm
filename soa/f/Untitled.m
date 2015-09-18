        % Channel response
        Gch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
            .*fiber.Hfiber(sim.f, tx);
        
        %'Fixed TD-FS-LE'
        Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); 

        if fiber.L*fiber.D(tx.lamb) ~= 0
            Hfiber = fiber.Hfiber(sim.f, tx)/fiber.link_attenuation(tx.lamb);
        else
            Hfiber = 1;
        end                

        Hmatched = Delay.*conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
            exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*Hfiber.*rx.elefilt.H(sim.f/sim.fs));
        % Note: Fiber frequency response is real, thus its group delay
        % is zero.
        
        
        % 'Adaptive TD-SR-LE'
                % Received Pulse Spectrum
        Hmatched = Delay.*conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
            exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*fiber.Hfiber(sim.f, tx));
        
        
                case 'Fixed TD-SR-LE'            
            Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); 

            if fiber.L*fiber.D(tx.lamb) ~= 0
                Hfiber = fiber.Hfiber(sim.f, tx)/fiber.link_attenuation(tx.lamb);
            else
                Hfiber = 1;
            end                
            
            Hmatched = Delay.*conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
                exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*Hfiber);
            % Note: Fiber frequency response is real, thus its group delay
            % is zero.