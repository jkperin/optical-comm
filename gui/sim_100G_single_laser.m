function sim_100G_single_laser
    close all, clc
    
    addpath f
    addpath ../f % general functions
    addpath ../mpam
    addpath ../soa
    addpath ../soa/f
    addpath ../apd
    addpath ../apd/f
    addpath ../ofdm
    addpath ../ofdm/f/

    getValue = @(h) str2double(get(h, 'String'));
    getLogicalValue = @(h) logical(get(h, 'Value'));
    
   %% GUI layout
   close all
   %  Create and then hide the GUI as it is being constructed.
   fSize_norm = [1 1];    
   f = figure('Visible','off', 'Units', 'Normalized', 'Position',[0 0 fSize_norm]);
 
   fSize = getpixelposition(f);
   fSize = fSize(3:4);
   
   panelWidht = 0.24;
   
   if fSize(2) < 1080
        HeaderFontSize = 10;
       FontSize = 9;
       BlockHeigth = 20/fSize(2);
   else
       HeaderFontSize = 12;
       FontSize = 11;
       BlockHeigth = 25/fSize(2);
   end
   
   %% Simulation Panel
   h.panel.sim_height = 0.06;
   maxY = 1 - h.panel.sim_height;
   h.panel.sim = uipanel('Title', 'Simulation', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY 1 h.panel.sim_height]);
   h.text.modulation = uicontrol('Parent', h.panel.sim, 'Style', 'text', 'String', 'Modulation Type:',...
       'Units', 'normalized', 'Position', [0 0 0.07 1], 'FontSize', FontSize);
   h.popup.modulation = uicontrol('Parent', h.panel.sim, 'Style','popupmenu', 'Units', 'normalized',...
       'String', {'M-PAM', 'M-CAP', 'DMT/OFDM'}, 'Position', [1 0 0.06 1],...
       'FontSize', FontSize-1, 'Callback',@popup_modulation_Callback);
   h.text.system = uicontrol('Parent', h.panel.sim, 'Style', 'text', 'String', 'System:',...
       'Units', 'normalized', 'Position', [1 0 0.05 1], 'FontSize', FontSize);
   h.popup.system = uicontrol('Parent', h.panel.sim, 'Style','popupmenu','Units', 'normalized',...
       'String', {'Basic', 'APD', 'SOA'}, 'Position', [1 0 0.06 1],...
       'FontSize', FontSize, 'Callback',@popup_system_Callback);
   h.text.results = uicontrol('Parent', h.panel.sim, 'Style', 'text', 'String', 'Results:',...
       'Units', 'normalized', 'Position', [1 0 0.05 1], 'FontSize', FontSize);
   h.popup.results = uicontrol('Parent', h.panel.sim, 'Style','popupmenu', 'Units', 'normalized',...
       'String', {'BER vs Transmitted Power', 'Power Penalty vs Modulator Cutoff Frequency',...
       'Power Penalty vs Fiber Length'}, 'Position', [1 0 0.15 1],...
       'FontSize', FontSize, 'Callback', @popup_results_Callback);
   h.button.run = uicontrol('Parent', h.panel.sim, 'Style', 'pushbutton', 'Units', 'normalized',...
       'String', 'Run Simulation', 'Position', [1-0.1 0 0.1 1],...
       'FontSize', FontSize, 'Callback', @run_Callback);
   align([h.text.modulation, h.popup.modulation, h.popup.system, h.text.results,...
       h.popup.results, h.text.system, h.button.run], 'None', 'Top');
   align([h.text.modulation, h.popup.modulation], 'Fixed', 10, 'None')
   align([h.popup.modulation, h.text.system], 'Fixed', 200, 'None')
   align([h.text.system, h.popup.system], 'Fixed', 10, 'None')
   align([h.popup.system, h.text.results], 'Fixed', 200, 'None')
   align([h.text.results, h.popup.results], 'Fixed', 10, 'None')
   

   %% Block diagram  
   h.panel.blockdiagram = uipanel('Title', 'Block Diagram', 'BackGroundColor', 'w', 'FontSize', HeaderFontSize,...
       'Position', [0.5 maxY-0.2 0.5 0.2]);
   h.axes.block_diagram = axes('Parent', h.panel.blockdiagram, 'Units', 'normalized',...
       'Position', [0 0.05 1 0.95], 'Color', 'w', 'Box', 'off');
   
   %% Results
   h.panel.results = uipanel('Title', 'Results', 'BackGroundColor', 'w', 'FontSize', HeaderFontSize,...
       'Position', [0.5 0 0.5 maxY-0.2]);
   h.axes.results = axes('Parent', h.panel.results, 'Units', 'Normalized', 'Position', [0.15 0.15 0.8 0.75], 'Box', 'on');
   h.clear = uicontrol('Parent', h.panel.results, 'Style', 'pushbutton', 'Units', 'normalized',...
       'String', 'Clear Plot', 'Position', [1-0.1 0.05 0.1 0.05],...
       'FontSize', FontSize, 'Callback', @clear_Callback);
   
   %% Transmitter Panel
   scale = 12.5; %9.5
   dH = 1/9;
   maxY = 1-h.panel.sim_height-scale*BlockHeigth;
   h.panel.tx = uipanel('Title', 'Transmitter', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY panelWidht scale*BlockHeigth]); 
   
   [h.text.Ptx, h.Ptx] = table_entry(h.panel.tx, [0 8*dH 1/scale], 'Transmitted Power Range (dBm):', '-30:2:-10'); 
   [h.text.fc, h.fc] = table_entry(h.panel.tx, [0 7*dH 1/scale], 'Modulator Cutoff Frequency (GHz):', 200);
   [h.text.lamb, h.lamb] = table_entry(h.panel.tx, [0 6*dH 1/scale], 'Wavelength (nm):', 1310);
   set(h.lamb, 'Callback', @(src, evt) update_dispersion_Callback(src, evt));
   
   [h.check.rex, h.text.rex, h.rex] = table_entry(h.panel.tx, [0 5*dH 1/scale], 'Extinction Ratio (dB):', -10, true,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
   
   [h.text.rin_variation, h.rin_variation] = table_entry(h.panel.tx, [0 2*dH 1/scale], 'RINmax - RINmin (dB/Hz)', 30, false);
   [h.text.rin_bw, h.rin_bw] = table_entry(h.panel.tx, [0 dH 1/scale], 'Noise Bandwidth for RIN Calculation (GHz)', 50, false);
   [h.check.rin_shape, h.text.rin_shape, h.rin_shape] = table_entry(h.panel.tx, [0 3*dH 1/scale], 'RIN Spectrum Shape Parameter (GHz)', 200, false,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, [internal_handles,...
       h.text.rin_variation, h.rin_variation h.text.rin_bw, h.rin_bw], false));
   [h.check.rin, h.text.rin, h.rin] = table_entry(h.panel.tx, [0 4*dH 1/scale], 'RIN (dB/Hz):', -140, true,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, [h.check.rin_shape internal_handles], false));
   
   [h.check.chirp, h.text.chirp, h.chirp] = table_entry(h.panel.tx, [0 dH/4 1/scale], 'Modulator chirp (alpha):', 2, false,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
   
   align([h.text.Ptx h.text.fc h.text.lamb h.check.rex h.check.rin  h.check.rin_shape h.text.rin_variation h.text.rin_bw  h.check.chirp], 'None', 'Fixed', 6);
   align([h.text.Ptx h.text.fc h.text.lamb h.text.rex h.text.rin h.text.rin_shape h.text.rin_variation  h.text.rin_bw h.text.chirp], 'Left', 'Fixed', 6)
   align([h.Ptx h.fc h.lamb h.rex h.rin h.rin_shape h.rin_variation h.rin_bw h.chirp], 'Right', 'Fixed', 6);
    
   %% Modulation Panel
   scale = 11.6;
   dH = 1/8;
   h.panel.mod = uipanel('Title', 'Modulation', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY-scale*BlockHeigth panelWidht scale*BlockHeigth]); 
   maxY = maxY - scale*BlockHeigth;
   
   [h.text.M, h.M] = table_entry(h.panel.mod, [0 7*dH 1/scale], 'Constellation Size:', 4);
   [h.text.pshape, h.pshape] = table_entry(h.panel.mod, [0 6*dH 1/scale], 'Pulse Shape:', 16);
   set(h.pshape, 'Style', 'popupmenu');
   set(h.pshape, 'String', {'Rectangular', 'Root-Raised Cosine'});
   set(h.pshape, 'Position', get(h.pshape, 'Position') + [-0.1 0 0.1 0])
   set(h.pshape, 'FontSize', FontSize-1);
   align([h.text.pshape, h.pshape], 'None', 'Center');
   
   [h.text.level, h.level] = table_entry(h.panel.mod, [0 5*dH 1/scale], 'Level Spacing:', 16);
   set(h.level, 'Style', 'popupmenu');
   set(h.level, 'String', {'Equal', 'Optimized'});
   set(h.level, 'Position', get(h.level, 'Position') + [-0.1 0 0.1 0])
   set(h.level, 'FontSize', FontSize-1);
   align([h.text.level, h.level], 'None', 'Center');
   
   [h.text.ros, h.ros] = table_entry(h.panel.mod, [0 4*dH 1/scale], 'Oversampling Ratio (ros):', 1);
   
   [h.text.Nc, h.Nc] = table_entry(h.panel.mod, [0 3*dH 1/scale], 'Number of Subcarriers: ', 64);
   [h.text.Nu, h.Nu] = table_entry(h.panel.mod, [0 2*dH 1/scale], 'Number of Used Subcarriers: ', 52);
   [h.text.palloc, h.palloc] = table_entry(h.panel.mod, [0 dH 1/scale], 'Power Allocation: ', 0);
   set(h.palloc, 'Style', 'popupmenu');
   set(h.palloc, 'String', {'Preemphasis', 'Opt. Bit Loading & Power allocation'});
   set(h.palloc, 'Position', get(h.palloc, 'Position') + [-0.1 0 0.1 0])
   set(h.palloc, 'FontSize', FontSize-1)
   align([h.text.palloc, h.palloc], 'None', 'Center');
   
   [h.text.rclip, h.rclip] = table_entry(h.panel.mod, [0 dH/4 1/scale], 'Clipping ratio (dB): ', 12);
      
   set_property([h.Nc, h.text.Nc, h.Nu, h.text.Nu, h.text.palloc, h.palloc, h.text.rclip, h.rclip], 'Enable', 'off');
   align([h.text.M h.text.pshape, h.text.level h.text.ros h.text.Nc h.text.Nu h.text.palloc, h.text.rclip], 'None', 'Fixed', 6);
   align([h.M h.pshape h.level h.ros h.Nc h.Nu h.palloc, h.rclip], 'Right', 'Fixed', 6)
   
   %% Fiber
   scale = 4.5; % 4
   dH = 1/3;
   h.panel.fiber = uipanel('Title', 'Fiber', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY-scale*BlockHeigth panelWidht scale*BlockHeigth]); 
   maxY = maxY - scale*BlockHeigth;
  
   [h.text.L, h.L] = table_entry(h.panel.fiber, [0 2*dH 1/scale], 'Length (km):', 0);
   [h.text.att, h.att] = table_entry(h.panel.fiber, [0 dH 1/scale], 'Attenuation (dB/km): ', 0.35);
   [h.text.D, h.D] = table_entry(h.panel.fiber, [0 dH/4 1/scale], 'Dispersion: (ps/(nm.km))', 0);
    align([h.text.L h.text.att h.text.D], 'None', 'Fixed', 6);
   align([h.L h.att h.D], 'Left', 'Fixed', 6)
   
   %% Receiver
   scale = 8.9; % 8
   dH = 1/6;
   maxY = maxY-scale*BlockHeigth;
   h.panel.rx = uipanel('Title', 'Receiver', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY panelWidht scale*BlockHeigth]); 
   
   [h.text.R, h.R] = table_entry(h.panel.rx, [0 6*dH 1/scale], 'Reponsivity (A/W): ', 1);
   [h.text.Id, h.Id] = table_entry(h.panel.rx, [0 5*dH 1/scale], 'Dark Current (nA):', 10);
   [h.check.shot, h.text.shot, temp1] = table_entry(h.panel.rx, [0 4*dH 1/scale], 'Shot Noise:', 0, true);
   set(temp1, 'Visible', 'off')
   [h.text.NEP, h.NEP] = table_entry(h.panel.rx, [0 3*dH 1/scale], 'NEP (pA/(Hz)^0.5):', 30);
   [h.text.rxfilter, h.rxfilterType] = table_entry(h.panel.rx, [0 2*dH 1/scale], 'Receiver Filter Order and Type', 0);
   [temp2, h.rxfilterN] = table_entry(h.panel.rx, [0 2*dH 1/scale], '', 4, false);
   set(temp2, 'Visible', 'off');
   [h.text.rxfilterBw, h.rxfilterBw] = table_entry(h.panel.rx, [0 dH/4 1/scale], 'Receiver Filter Bandwidth (GHz):', 30, false);
  
   align([h.text.R, h.text.Id, h.text.shot, h.text.NEP, h.text.rxfilter, h.text.rxfilterBw], 'None', 'Fixed', 6);
   align([h.R, h.Id, temp1, h.NEP,  h.rxfilterN, h.rxfilterBw], 'None', 'Fixed', 6);
   align([h.text.rxfilter, h.rxfilterType, h.rxfilterN], 'None', 'Bottom')
  
   set(h.rxfilterN, 'Position', get(h.rxfilterN, 'Position') + [-0.26 0 -0.05 0])
   set(h.rxfilterType, 'Position', get(h.rxfilterType, 'Position') + [-0.1 0 0.1 0])
   set(h.rxfilterType, 'Style', 'popupmenu');
   set(h.rxfilterType, 'String', {'Matched', 'Gaussian', 'Bessel'});
   set(h.rxfilterType, 'FontSize', FontSize);
   
   set(h.rxfilterType, 'Callback', @(src,evt) disable_handles_Callback(src, evt, [h.rxfilterN  h.rxfilterBw h.text.rxfilterBw], 'Matched'))
   
   
   align([h.text.rxfilter, h.rxfilterType], 'None', 'Center');
   align([h.text.shot, h.check.shot], 'None', 'Bottom')
   
   %% APD
   scale = 4.5;
   dH = 1/3;
   maxY = 1-h.panel.sim_height-scale*BlockHeigth;
   h.panel.apd = uipanel('Title', 'APD', 'FontSize', HeaderFontSize,...
       'Position', [panelWidht+0.01 maxY panelWidht scale*BlockHeigth]); 
   
   [h.check.GBw, h.text.GBw, h.GBw] = table_entry(h.panel.apd, [0 2*dH 1/scale], 'Finite Gain-Bandwidth Product (GHz):', 340, true,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
   [h.text.ka, h.ka] = table_entry(h.panel.apd, [0 dH 1/scale], 'Impact Ionization Factor (ka):', 0.09);
   [h.check.Gapd, h.text.Gapd, h.Gapd] = table_entry(h.panel.apd, [0 dH/4 1/scale], 'Set Gain (dB):', 10, false,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
   
   align([h.check.GBw, h.text.ka, h.check.Gapd], 'None', 'Fixed', 6);
   align([h.text.GBw, h.text.ka, h.text.Gapd], 'None', 'Fixed', 6);
   align([h.GBw, h.ka, h.Gapd], 'None', 'Fixed', 6);   
   
   %% SOA
   scale = 7;
   dH = 1/5;
   maxY = maxY - scale*BlockHeigth;
   h.panel.soa = uipanel('Title', 'SOA', 'FontSize', HeaderFontSize,...
       'Position', [panelWidht+0.01 maxY panelWidht scale*BlockHeigth]); 
   
   [h.text.Gsoa, h.Gsoa] = table_entry(h.panel.soa, [0 4*dH 1/scale], 'Gain (dB):', 20);
   [h.text.Fn, h.Fn] = table_entry(h.panel.soa, [0 3*dH 1/scale], 'Noise Figure (dB):', 7);
      
   [h.text.optfiltType, h.optfiltType] = table_entry(h.panel.soa, [0 2*dH 1/scale], 'Optical Filter Type:', 0);
    set(h.optfiltType, 'Style', 'popupmenu');
    set(h.optfiltType, 'String', {'Fiber Brag Gratting', 'Fabry-Perot'});
    set(h.optfiltType, 'Position', get(h.optfiltType, 'Position') + [-0.1 0 0.1 0])
    set(h.optfiltType, 'FontSize', FontSize-1)
    align([h.text.optfiltType, h.optfiltType], 'None', 'Center');
    
    [h.text.Bopt, h.Bopt] = table_entry(h.panel.soa, [0 dH 1/scale], 'Optical Filter Bandwidth (GHz):', 200);
    [h.check.polarizer, h.text.polarizer, temp] = table_entry(h.panel.soa, [0 dH/4 1/scale], 'Polarization Tracking', 0, true,...
        @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
    set(temp, 'Visible', 'off')
     
   align([h.text.Gsoa, h.text.Fn, h.text.optfiltType, h.text.Bopt, h.text.polarizer], 'None', 'Fixed', 6);
   align([h.Gsoa, h.Fn, h.optfiltType, h.Bopt, h.check.polarizer], 'None', 'Fixed', 6); 
         
   %% Simulation
   scale = 9;
   dH = 1/6;
   maxY = maxY - scale*BlockHeigth;
   h.panel.param = uipanel('Title', 'Simulation', 'FontSize', HeaderFontSize,...
       'Position', [panelWidht+0.01 maxY panelWidht scale*BlockHeigth]); 
   
   [h.text.Rb, h.Rb] = table_entry(h.panel.param, [0 5*dH 1/scale], 'Bit Rate (Gbps):', 100);
   [h.text.BER, h.BER] = table_entry(h.panel.param, [0 4*dH 1/scale], 'Target BER: ', 1e-4);
   [h.text.Nsymb, h.Nsymb] = table_entry(h.panel.param, [0 3*dH 1/scale], 'Number of Symbols:', 2^15);
   [h.text.Mct, h.Mct] = table_entry(h.panel.param, [0 2*dH 1/scale], 'Oversampling Ratio for Continuous Time:', 15);
   
   [h.check.quantiz, h.text.ENOB, h.ENOB] = table_entry(h.panel.param, [0 dH 1/scale], 'ENOB (bits):', 6, false,...
       @(src, evt, internal_handles) disable_handles_Callback(src, evt, internal_handles, false));
      
   [h.text.Lseq, h.Lseq] = table_entry(h.panel.param, [0 dH/4 1/scale], 'de Bruijn Sub-Sequence Length:', 2);
   align([h.text.Rb, h.text.BER, h.text.Nsymb, h.text.Mct, h.text.ENOB, h.text.Lseq], 'None', 'Fixed', 6);
   align([h.text.Rb, h.text.BER, h.text.Nsymb, h.text.Mct, h.check.quantiz, h.text.Lseq], 'None', 'Fixed', 6);
   align([h.Rb, h.BER, h.Nsymb, h.Mct, h.ENOB h.Lseq], 'None', 'Fixed', 6);
   
   %Create a plot in the axes.
   popup_system_Callback(h.popup.system, 0)
%    current_data = peaks_data;
%    surf(current_data);
   % Assign the GUI a name to appear in the window title.
   set(f, 'Name', '100G single-laser simulation');
   % Move the GUI to the center of the screen.
   movegui(f,'center')
   % Make the GUI visible.
   set(f, 'Visible', 'on');
 
   %% Pushbutton Callbacks.
    function clear_Callback(source, eventdata)
        cla(h.axes.results)
    end
        
    function run_Callback(source, eventdata)
        clc
        
        [mpam, ofdm1, tx, fiber1, soa1, apd1, rx, sim] = build_simulation(h);
        
        switch getOption(h.popup.results)
            case 'BER vs Transmitted Power'
                switch getOption(h.popup.modulation)
                    case 'M-PAM'
                        [ber, GdB] = mpam_ber_vs_Ptx(mpam, tx, fiber1, soa1, apd1, rx, sim);
                    case 'DMT/OFDM'
                        ber = ofdm_ber_vs_Ptx(ofdm1, tx, fiber1, rx, sim);
                end
                
                if isempty(soa1)
                    handle = h.Gapd;
                else
                    handle = h.Gsoa;
                end
                
                if exist('GdB', 'var') && ~isempty(GdB)
                    set(handle, 'String', num2str(GdB, 2));
                end

                axes(h.axes.results), hold on, box on, grid on
                hline = plot(tx.PtxdBm, log10(ber.count), '--o');
                plot(tx.PtxdBm, log10(ber.est), '-', 'Color', get(hline, 'Color'))
                xlabel('Transmitted Power (dBm)')
                ylabel('log_{10}(BER)')
                legend('Montecarlo', 'Estimated')
                a = axis;
                axis([min(a(1), tx.PtxdBm(1)) max(a(2), tx.PtxdBm(end)) ceil(2*log10(sim.BERtarget)) 0])
                
            case 'Power Penalty vs Modulator Cutoff Frequency'
                switch getOption(h.popup.modulation)
                    case 'M-PAM'
                        %pp = mpam_pp_vs_mod_cutoff(mpam, tx, fiber1, soa1, apd1, rx, sim);
                        error('Not implemented yet.')
                    case 'DMT/OFDM'
                        pp = ofdm_pp_vs_mod_cutoff(ofdm1, tx, fiber1, rx, sim);
                end
                
                axes(h.axes.results), hold on, box on, grid on
                hline = plot(tx.modulator.Fc/1e9, pp.power_pen_ook_m, '--o', 'LineWidth', 2);
                plot(tx.modulator.Fc/1e9, pp.power_pen_ook_e, '-', 'LineWidth', 2, 'Color', get(hline, 'Color'))
                legend('Montecarlo', 'Estimated', 'Location', 'NorthEast')
                xlabel('Modulator Cutoff Frequency (GHz)')
                ylabel('Power Penalty w.r. to NRZ-OOK (dB)')
                
            case 'Power Penalty vs Fiber Length'
                switch getOption(h.popup.modulation)
                    case 'M-PAM'
%                         pp = mpam_pp_vs_L(mpam, tx, fiber1, soa1, apd1, rx, sim);
                          error('Not implemented yet.')
                    case 'DMT/OFDM'
                        pp = ofdm_pp_vs_L(ofdm1, tx, fiber1, rx, sim);
                end
                
                axes(h.axes.results), hold on, box on, grid on
                hline = plot(sim.fiberL/1e3, pp.power_pen_ook_m, '--o', 'LineWidth', 2);
                plot(sim.fiberL/1e3, pp.power_pen_ook_e, '-', 'LineWidth', 2, 'Color', get(hline, 'Color'))
                legend('Montecarlo', 'Estimated', 'Location', 'NorthEast')
                xlabel('Fiber Length (km)')
                ylabel('Power Penalty w.r. to NRZ-OOK (dB)')
                
            otherwise
                error('Selected results not implemented yet')
        end             
    end   
   
   %% Popup Callbacks
    function popup_modulation_Callback(source, eventdata) 
         % Determine the selected data set.
         str = get(source, 'String');
         val = get(source, 'Value');
         % Set current data to the selected data set.
         if strcmp(str{val}, 'M-PAM') % User selects Sinc.
             set_property([h.Nc, h.level, h.text.level, h.text.Nc, h.Nu,...
                 h.text.Nu, h.palloc, h.text.palloc, h.text.rclip, h.rclip], 'Enable', 'off')
             set_property([h.level, h.text.level h.pshape h.text.pshape, h.text.rxfilter, h.rxfilterN,...
                 h.text.rxfilterBw, h.rxfilterBw, h.text.ros, h.ros, h.Lseq, h.text.Lseq], 'Enable', 'on')
                          
             set(h.popup.system, 'String', {'Basic', 'SOA', 'APD'})
             set(h.rxfilterType, 'String', {'Matched', 'Gaussian', 'Bessel'})
             set(h.rxfilterType, 'Value', 1);
             set(h.rxfilterType, 'Callback', @(src,evt) disable_handles_Callback(src, evt, [h.rxfilterN  h.rxfilterBw h.text.rxfilterBw], 'Matched'))
             h.rxfilterType.Callback(h.rxfilterType, 0);
             
             switch getOption(h.popup.system)
                 case 'Basic'
                     imshow('data/mpam-basic.png', 'Parent', h.axes.block_diagram)
                 case 'SOA'
                     imshow('data/mpam-soa.png', 'Parent', h.axes.block_diagram)
                 case 'APD'
                     imshow('data/mpam-apd.png', 'Parent', h.axes.block_diagram)
             end
             
             set(h.M, 'String', 4); 
             set(h.ros, 'String', 1);
             set(h.Ptx, 'String', '-30:2:-10')
             set(h.Nsymb, 'String', num2str(2^14))
             set(h.Mct, 'String', num2str(15))     
         else
             set_property([h.Nc, h.text.Nc, h.Nu, h.text.Nu, h.palloc, h.text.palloc,...
                 h.text.rclip, h.rclip, h.text.rxfilter, h.rxfilterN,...
                 h.text.rxfilterBw, h.rxfilterBw], 'Enable', 'on')
             set_property([h.level, h.text.level h.pshape h.text.pshape,...
                 h.text.ros, h.ros, h.Lseq, h.text.Lseq, h.text.rxfilterBw, h.rxfilterBw], 'Enable', 'off')
             
             set(h.rxfilterType, 'String', {'Gaussian', 'Bessel'})
             set(h.rxfilterType, 'Value', 1);
             set(h.rxfilterType, 'Callback', @(src,evt) [])
             
             set(h.popup.system, 'String', {'Basic'})
            
             set(h.ros, 'String', num2str(getValue(h.Nc)/getValue(h.Nu), 3))
                      
             imshow('data/ofdm.png', 'Parent', h.axes.block_diagram)
             
             set(h.M, 'String', 16); 
             set(h.Ptx, 'String', '-20:2:-4')
             set(h.Nsymb, 'String', num2str(2^12))
             set(h.Mct, 'String', num2str(9))             
         end
    end

    function popup_results_Callback(source, eventdata) 
         str = get(source, 'String');
         val = get(source, 'Value');
         % Set current data to the selected data set.
         switch str{val}
             case 'BER vs Transmitted Power'
                 set_property([h.text.Ptx h.Ptx], 'Enable', 'on');
                 set(h.text.Ptx, 'String', 'Transmitted Power Range (dBm):')
                 set(h.Ptx, 'String', '-20:2:-4')
                 
                 set(h.text.fc, 'String', 'Modulator Cutoff Frequency (GHz):')
                 set(h.fc, 'String', '30')
                 
                 set(h.text.L, 'String', 'Fiber Length (km):')
                 set(h.L, 'String', '0')
                                 
             case 'Power Penalty vs Modulator Cutoff Frequency'
                 set_property([h.text.Ptx h.Ptx], 'Enable', 'off');
                 
                 set(h.text.fc, 'String', 'Modulator Cutoff Frequency Range (GHz):')
                 set(h.fc, 'String', '20:5:50')
                 
                 set(h.text.L, 'String', 'Fiber Length (km):')
                 set(h.L, 'String', '0')
                 
             case 'Power Penalty vs Fiber Length'
                 set_property([h.text.Ptx h.Ptx], 'Enable', 'off');
                 
                 set(h.text.fc, 'String', 'Modulator Cutoff Frequency (GHz):')
                 set(h.fc, 'String', '30')
                 
                 set(h.text.L, 'String', 'Fiber Length Range (km):')
                 set(h.L, 'String', '0:5')
                 
                 set(h.check.chirp, 'Value', true)
                 h.check.chirp.Callback(h.check.chirp, 0)
         end               
    end

    function update_dispersion_Callback(src, evt)       
        lamb = 1e-9*str2double(get(h.lamb, 'String'));
        
        D = 1e6*fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); % Dispersion curve
        
        set(h.D, 'String', num2str(D))        
    end

    % Disable handles if get(source, 'String') == value
    function disable_handles_Callback(source, ~, handles, value)
         str = get(source, 'String');
         val = get(source, 'Value');
         
         off = false;
         if islogical(value)
             off = logical(val) == value;
         else
             off = strcmp(str{val}, value);
         end
             
         if off
             s = 'off';
         else
             s = 'on';
         end
         
          for k = 1:length(handles)
             set(handles, 'Enable', s)
         end
    end
        
    function popup_system_Callback(source, ~)
         % Determine the selected data set.
         str = get(source, 'String');
         val = get(source, 'Value');
         % Set current data to the selected data set.
         switch str{val}
             case 'Basic'
                 apd_e = 'off';
                 soa_e = 'off';
                 
                 set(h.popup.modulation, 'String', {'M-PAM', 'DMT/OFDM'})
                 
                 if strcmp(getOption(h.popup.modulation), 'M-PAM')
                    imshow('data/mpam-basic.png', 'Parent', h.axes.block_diagram)
                 else
                     imshow('data/ofdm.png', 'Parent', h.axes.block_diagram)
                 end
                 
             case 'APD'
                 apd_e = 'on';
                 soa_e = 'off';
                 
                 set(h.popup.modulation, 'String', {'M-PAM'})
                                  
                 imshow('data/mpam-apd.png', 'Parent', h.axes.block_diagram)
             case 'SOA'
                 apd_e = 'off';
                 soa_e = 'on';
                 
                 set(h.popup.modulation, 'String', {'M-PAM'})
                 
                 imshow('data/mpam-soa.png', 'Parent', h.axes.block_diagram)
         end
         
         apd_child = allchild(h.panel.apd);
         soa_child = allchild(h.panel.soa);
         for k = 1:length(soa_child)
             set(soa_child(k), 'Enable', soa_e);
         end
         for k = 1:length(apd_child)
             set(apd_child(k), 'Enable', apd_e);
         end
         
         %% Call callbacks
         switch str{val}
             case 'APD'
                h.check.GBw.Callback(h.check.GBw, []);
                h.check.Gapd.Callback(h.check.Gapd, []);
             case 'SOA'
                h.check.polarizer.Callback(h.check.polarizer, []);
         end
             
    end

    %% Auxiliary functions
    function set_property(handle_array, property, value)
        for k = 1:length(handle_array)
            set(handle_array, property, value)
        end
    end
   
    function varargout = table_entry(parent, norm_pos, str, value, checked, Callback)
        if nargin < 5
            checked = true;
        end
        
        if nargin < 6
            Callback = [];
        end
        
        if length(norm_pos) == 2
            norm_pos(3) = 0.1;
        end
        
        cb_width = 0.05;
        edit_width = 0.2;
        
        checkbox = uicontrol('Parent', parent, 'Style', 'checkbox',...
            'String', '', 'Units', 'Normalized','FontSize', FontSize, 'Value', checked);
        cb_pos = get(checkbox, 'Position');
        set(checkbox, 'Position', [norm_pos(1:2) cb_width norm_pos(3)]);

        text = uicontrol('Parent', parent, 'Style', 'text', 'HorizontalALignment', 'left',...
        'String', str, 'Units', 'Normalized', 'FontSize', FontSize,...
        'Position', [norm_pos(1)+cb_width norm_pos(2) 1-(norm_pos(1)+cb_width+edit_width+0.13) norm_pos(3)]); 

        edit = uicontrol('Parent', parent, 'Style', 'edit', 'String', num2str(value),...
           'HorizontalALignment', 'right', 'Units', 'Normalized', 'FontSize', FontSize,...
           'Position', [1-edit_width-0.01, norm_pos(2), edit_width norm_pos(3)]); 

       if ~checked
           set(text, 'Enable', 'off')
           set(edit, 'Enable', 'off')
       end
       
        if ~isempty(Callback)
            set(checkbox, 'Callback', @(src, evt) Callback(src, evt, [text edit]))
        end
       
        uistack(edit, 'top')
        align([text, checkbox, edit], 'None', 'Center')
        if nargout == 3
            varargout = {checkbox, text, edit};
        else
            varargout = {text, edit};
            set(checkbox, 'Visible', 'off');
        end
    end           
end

function str = getOption(h)
    list = get(h, 'String');
    str = list{get(h, 'Value')};
end

