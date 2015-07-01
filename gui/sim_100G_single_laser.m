function sim_100G_single_laser2
    close all
   %  Create and then hide the GUI as it is being constructed.
   fSize_norm = [1 1];    
   f = figure('Visible','off', 'Units', 'Normalized', 'Position',[0 0 fSize_norm]);
 
   fSize = getpixelposition(f);
   fSize = fSize(3:4);
   
   BlockHeigth = 30/fSize(2);
   if fSize(2) <= 1080
        HeaderFontSize = 10;
       FontSize = 9;
   else
       HeaderFontSize = 12;
       FontSize = 11;
   end
   
   % Simulation Panel
   hpanel_sim_height = 0.07;
   maxY = 1 - hpanel_sim_height;
   hpanel_sim = uipanel('Title', 'Simulation', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY 1 hpanel_sim_height]);
   htext_modulation = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'Modulation Type:',...
       'Units', 'normalized', 'Position', [0 0 0.07 1], 'FontSize', FontSize);
   hpopup_modulation = uicontrol('Parent', hpanel_sim, 'Style','popupmenu', 'Units', 'normalized',...
       'String', {'M-PAM', 'M-CAP', 'DMT/OFDM'}, 'Position', [1 0 0.06 1],...
       'FontSize', FontSize-1, 'Callback',@popup_modulation_Callback);
   htext_system = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'System:',...
       'Units', 'normalized', 'Position', [1 0 0.05 1], 'FontSize', FontSize);
   hpopup_system = uicontrol('Parent', hpanel_sim, 'Style','popupmenu','Units', 'normalized',...
       'String', {'Basic', 'APD', 'SOA'}, 'Position', [1 0 0.06 1],...
       'FontSize', FontSize, 'Callback',@popup_system_Callback);
   htext_results = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'Results:',...
       'Units', 'normalized', 'Position', [1 0 0.05 1], 'FontSize', FontSize);
   hpopup_results = uicontrol('Parent', hpanel_sim, 'Style','popupmenu', 'Units', 'normalized',...
       'String', {'BER vs Transmitted Power', 'Power Penalty vs Modulator Cutoff Frequency',...
       'SOA/APD Gain'}, 'Position', [1 0 0.15 1],...
       'FontSize', FontSize, 'Callback', @popup_results_Callback);
   hrun = uicontrol('Parent', hpanel_sim, 'Style', 'pushbutton', 'Units', 'normalized',...
       'String', 'Run Simulation', 'Position', [1-0.1 0 0.1 1],...
       'FontSize', FontSize, 'Callback', @run_Callback);
   align([htext_modulation, hpopup_modulation, hpopup_system, htext_results,...
       hpopup_results, htext_system, hrun], 'None', 'Top');
   align([htext_modulation, hpopup_modulation], 'Fixed', 10, 'None')
   align([hpopup_modulation, htext_system], 'Fixed', 200, 'None')
   align([htext_system, hpopup_system], 'Fixed', 10, 'None')
   align([hpopup_system, htext_results], 'Fixed', 200, 'None')
   align([htext_results, hpopup_results], 'Fixed', 10, 'None')
   
   
   
   %% Block diagram  
   hpanel_blockdiagram = uipanel('Title', 'Block Diagram', 'FontSize', HeaderFontSize,...
       'Position', [0.5 1-hpanel_sim_height-.2 0.5 .2]);
   block_diagram = axes('Parent', hpanel_blockdiagram, 'Units', 'Pixels',...
       'Position', [5 5 fSize(1) 120], 'Color', 'w', 'Box', 'off');
   
   %% Results
   hpanel_results = uipanel('Title', 'Results', 'FontSize', HeaderFontSize,...
       'Position', [0.5 0 0.5 1-hpanel_sim_height-.2]);
   results = axes('Parent', hpanel_results, 'Units', 'Normalized', 'Position', [0.15 0.15 0.8 0.8], 'Box', 'on');
   
   %% Modulation Panel
   %% Transmitter Panel
   scale = 8;
   dH = 1/6;
   hpanel_tx = uipanel('Title', 'Transmitter', 'FontSize', HeaderFontSize,...
       'Position', [0 1-hpanel_sim_height-scale*BlockHeigth 0.35 scale*BlockHeigth]); 
   maxY = 1-hpanel_sim_height-scale*BlockHeigth;
   
   [htext_lambda, hlambda] = table_entry(hpanel_tx, [0 6*dH 1/scale], 'Wavelength (nm):', 1310);
   [htext_fc, hfc] = table_entry(hpanel_tx, [0 5*dH 1/scale], 'Modulator Cutoff Frequency (GHz):', 30);
   [htext_kappa, hkappa] = table_entry(hpanel_tx, [0 4*dH 1/scale], 'Insertion Loss (dB):', 0);
   [hcheck_rex, htext_rex, hrex] = table_entry(hpanel_tx, [0 3*dH 1/scale], 'Extinction Ratio (dB):', -15, true);
   [hcheck_rin, htext_rin, hrin] = table_entry(hpanel_tx, [0 2*dH 1/scale], 'RIN (dB/Hz):', -150, true);
   [hcheck_chirp, htext_chirp, hchirp] = table_entry(hpanel_tx, [0 dH/4 1/scale], 'Modulator chirp (alpha):', -2);
   align([htext_lambda htext_fc htext_kappa hcheck_rex hcheck_rin hcheck_chirp], 'None', 'Fixed', 6);
   align([htext_lambda htext_fc htext_kappa htext_rex htext_rin htext_chirp], 'Left', 'Fixed', 6)
   align([hlambda hfc hkappa hrex hrin hchirp], 'Right', 'Fixed', 6);
    
   %% Modulation Panel
   scale = 6;
   dH = 1/3;
   hpanel_mod = uipanel('Title', 'Modulation', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY-scale*BlockHeigth 0.35 scale*BlockHeigth]); 
   maxY = maxY - scale*BlockHeigth;
   
   [htext_M, hM] = table_entry(hpanel_mod, [0 2*dH 1/scale], 'Constellation Size:', 16);
   [htext_Nc, hNc] = table_entry(hpanel_mod, [0 dH 1/scale], 'Number of Subcarriers: ', 64);
   [htext_Nu, hNu] = table_entry(hpanel_mod, [0 dH/4 1/scale], 'Number of Used Subcarriers: ', 52);
   [htext_palloc, hpalloc] = table_entry(hpanel_mod, [0 dH/4 1/scale], 'Power Allocation: ', 0);
   set(hpalloc, 'Style', 'popupmenu');
   set(hpalloc, 'String', {'Preemphasis', 'Opt. Bit Loading & Power allocation'});
   set(hpalloc, 'Position', get(hpalloc, 'Position') + [-0.1 0 0.1 0])
   set(hpalloc, 'FontSize', FontSize-1)
   align([htext_palloc, hpalloc], 'None', 'Center');
   set_property([hNc, htext_Nc, hNu, htext_Nu, htext_palloc, hpalloc], 'Enable', 'off');
   align([htext_M htext_Nc htext_Nu htext_palloc], 'None', 'Fixed', 6);
   align([hM hNc hNu hpalloc], 'Right', 'Fixed', 6)
   
   %% Fiber
   scale = 4;
   dH = 1/3;
   hpanel_fiber = uipanel('Title', 'Fiber', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY-scale*BlockHeigth 0.35 scale*BlockHeigth]); 
   maxY = maxY - scale*BlockHeigth;
  
   [htext_L, hL] = table_entry(hpanel_fiber, [0 2*dH 1/scale], 'Length (km):', 0);
   [htext_att, hatt] = table_entry(hpanel_fiber, [0 dH 1/scale], 'Attenuation (dB/km): ', 0.35);
   [htext_D, hD] = table_entry(hpanel_fiber, [0 dH/4 1/scale], 'Dispersion: (ps/(nm.km)', 0);
    align([htext_L htext_att htext_D], 'None', 'Fixed', 6);
   align([hL hatt hD], 'Left', 'Fixed', 6)
   
   %% Receiver
   scale = 7;
   dH = 1/6;
   maxY = maxY-scale*BlockHeigth;
   hpanel_rx = uipanel('Title', 'Receiver', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY 0.35 scale*BlockHeigth]); 
   
   [htext_NEP, hNEP] = table_entry(hpanel_rx, [0 5*dH 1/scale], 'NEP (pA/(Hz)^0.5):', 30);
   [htext_R, hR] = table_entry(hpanel_rx, [0 4*dH 1/scale], 'Reponsivity (A/W): ', 1);
   [htext_Id, hId] = table_entry(hpanel_rx, [0 3*dH 1/scale], 'Dark Current (nA):', 10);
   [htext_rxfilter, hrxfilterType] = table_entry(hpanel_rx, [0 2*dH 1/scale], 'Receiver Filter Order and Type', 0);
   [temp, hrxfilterN] = table_entry(hpanel_rx, [0 2*dH 1/scale], '', 4);
   set(temp, 'Visible', 'off');
   [htext_rxfilterBw, hrxfilterBw] = table_entry(hpanel_rx, [0 dH/4 1/scale], 'Receiver Filter Bandwidth (GHz):', 30);
   align([htext_NEP, htext_R, htext_Id, htext_rxfilter, htext_rxfilterBw], 'None', 'Fixed', 6);
   align([hNEP, hR, hId, hrxfilterN, hrxfilterBw], 'None', 'Fixed', 6);
   align([htext_rxfilter, hrxfilterType, hrxfilterN], 'None', 'Bottom')
   set(hrxfilterN, 'Position', get(hrxfilterN, 'Position') + [-0.20 0 -0.05 0])
   set(hrxfilterType, 'Position', get(hrxfilterType, 'Position') + [-0.1 0 0.1 0])
   set(hrxfilterType, 'Style', 'popupmenu');
   set(hrxfilterType, 'String', {'Gaussian', 'Bessel', 'Matched'});
   set(hrxfilterType, 'FontSize', FontSize-1);
   align([htext_rxfilter, hrxfilterType], 'None', 'Center');
   
   %% Simulation
   scale = 7;
   dH = 1/6;
   maxY = maxY-scale*BlockHeigth;
   hpanel_param = uipanel('Title', 'Receiver', 'FontSize', HeaderFontSize,...
       'Position', [0 maxY 0.35 scale*BlockHeigth]); 
   
   [htext_Rb, hRb] = table_entry(hpanel_param, [0 4*dH 1/scale], 'Bit Rate (Gbps):', 100);
   [htext_BER, hBER] = table_entry(hpanel_param, [0 3*dH 1/scale], 'Target BER: ', 1e-4);
   [htext_Nsymb, hNsymb] = table_entry(hpanel_param, [0 2*dH 1/scale], 'Number of Symbols:', 2^15);
   [htext_Mct, hMct] = table_entry(hpanel_param, [0 dH 1/scale], 'Oversampling Ratio for Continuous Time:', 8);
   [hcheck_quantiz, htext_ENOB, hENOB] = table_entry(hpanel_param, [0 dH/4 1/scale], 'ENOB (bits):', 6);
   align([htext_Rb, htext_BER, htext_Nsymb, htext_Mct, htext_ENOB], 'None', 'Fixed', 6);
   align([hRb, hBER, hNsymb, hMct, hENOB], 'None', 'Fixed', 6);
  
   
   
%    %  Construct the components.
%    hsurf = uicontrol('Style','pushbutton','String','Surf',...
%           'Position',[315,220,70,25],...
%           'Callback',@surfbutton_Callback);
%    hmesh = uicontrol('Style','pushbutton','String','Mesh',...
%           'Position',[315,180,70,25],...
%           'Callback',@meshbutton_Callback);
%    hcontour = uicontrol('Style','pushbutton',...
%           'String','Countour',...
%           'Position',[315,135,70,25],...
%           'Callback',@contourbutton_Callback); 
%    htext = uicontrol('Style','text','String','Select Data',...
%           'Position',[325,90,60,15]);
%    hpopup = uicontrol('Style','popupmenu',...
%           'String',{'Peaks','Membrane','Sinc'},...
%           'Position',[300,50,100,25],...
%           'Callback',@popup_menu_Callback);
%    ha = axes('Units','Pixels','Position',[50,60,200,185]); 
%    align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
%    
   % Create the data to plot.
   peaks_data = peaks(35);
   membrane_data = membrane;
   [x,y] = meshgrid(-8:.5:8);
   r = sqrt(x.^2+y.^2) + eps;
   sinc_data = sin(r)./r;
   
   % Initialize the GUI.
   % Change units to normalized so components resize 
   % automatically.
   f.Units = 'normalized';
   ha.Units = 'normalized';
   hsurf.Units = 'normalized';
   hmesh.Units = 'normalized';
   hcontour.Units = 'normalized';
   htext.Units = 'normalized';
   hpopup.Units = 'normalized';
   
   %Create a plot in the axes.
%    current_data = peaks_data;
%    surf(current_data);
   % Assign the GUI a name to appear in the window title.
   f.Name = 'Simple GUI';
   % Move the GUI to the center of the screen.
   movegui(f,'center')
   % Make the GUI visible.
   f.Visible = 'on';
 
    function popup_modulation_Callback(source, eventdata) 
         % Determine the selected data set.
         str = source.String;
         val = source.Value;
         % Set current data to the selected data set.
         if ~strcmp(str{val}, 'DMT/OFDM') % User selects Sinc.
             set_property([hNc, htext_Nc, hNu, htext_Nu, hpalloc, htext_palloc], 'Enable', 'off')
         else
             set_property([hNc, htext_Nc, hNu, htext_Nu, hpalloc, htext_palloc], 'Enable', 'on')
         end
    end
  

    function set_property(handle_array, property, value)
        for k = 1:length(handle_array)
            set(handle_array, property, value)
        end
    end
   
    function varargout = table_entry(parent, norm_pos, str, value, checked)
        if nargin < 5
            checked = false;
        end
        
        if length(norm_pos) == 2
            norm_pos(3) = 0.1;
        end
        
      checkbox = uicontrol('Parent', parent, 'Style', 'checkbox',...
       'String', '', 'Units', 'Normalized','FontSize', FontSize, 'Value', checked);
       cb_pos = get(checkbox, 'Position');
       set(checkbox, 'Position', [norm_pos(1:2) 0.03 norm_pos(3)]);
      text = uicontrol('Parent', parent, 'Style', 'text', 'HorizontalALignment', 'left',...
       'String', str, 'Units', 'Normalized', 'Position', [norm_pos(1)+0.03 norm_pos(2) 0.8 norm_pos(3)],...
       'FontSize', FontSize); % wavelength
      edit = uicontrol('Parent', parent, 'Style', 'edit', 'String', num2str(value),...
           'HorizontalALignment', 'right', 'Units', 'Normalized',...
           'Position', [norm_pos(1)+0.85, norm_pos(2), 0.14 norm_pos(3)],...
           'FontSize', FontSize); % wavelength
      align([text, checkbox, edit], 'None', 'Center')
        if nargout == 3
            varargout = {checkbox, text, edit};
        else
            varargout = {text, edit};
            set(checkbox, 'Visible', 'off');
        end
    end
   
   %  Callbacks for simple_gui. These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 
   %  Pop-up menu callback. Read the pop-up menu Value property
   %  to determine which item is currently displayed and make it
   %  the current data.
      function popup_menu_Callback(source,eventdata) 
         % Determine the selected data set.
         str = source.String;
         val = source.Value;
         % Set current data to the selected data set.
         switch str{val};
         case 'Peaks' % User selects Peaks.
            current_data = peaks_data;
         case 'Membrane' % User selects Membrane.
            current_data = membrane_data;
         case 'Sinc' % User selects Sinc.
            current_data = sinc_data;
         end
      end
  
   % Push button callbacks. Each callback plots current_data in
   % the specified plot type.
 
   function surfbutton_Callback(source,eventdata) 
   % Display surf plot of the currently selected data.
      surf(current_data);
   end
 
   function meshbutton_Callback(source,eventdata) 
   % Display mesh plot of the currently selected data.
      mesh(current_data);
   end
 
   function contourbutton_Callback(source,eventdata) 
   % Display contour plot of the currently selected data.
      contour(current_data);
   end 
 
end 