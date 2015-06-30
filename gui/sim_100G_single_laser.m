function sim_100G_single_laser2
    close all
   %  Create and then hide the GUI as it is being constructed.
   fSize = [1000 700];
   f = figure('Visible','off','Position',[360,500, fSize]);
 
   HeaderFontSize = 12;
   FontSize = 11;
   
   % Simulation Panel
   hpanel_sim_height = 0.08;
   hpanel_sim = uipanel('Title', 'Simulation', 'FontSize', HeaderFontSize,...
       'Position', [0 1-hpanel_sim_height 1 hpanel_sim_height]);
   htext_modulation = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'Modulation Type:',...
       'Position', [0 5 120 20], 'FontSize', FontSize);
   hpopup_modulation = uicontrol('Parent', hpanel_sim, 'Style','popupmenu',...
       'String', {'M-PAM', 'M-CAP', 'DMT/OFDM'}, 'Position', [125 9 100 20],...
       'FontSize', FontSize, 'Callback',@popup_modulation_Callback);
   htext_system = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'System:',...
       'Position', [400 5 60 20], 'FontSize', FontSize);
   hpopup_system = uicontrol('Parent', hpanel_sim, 'Style','popupmenu',...
       'String', {'Basic', 'APD', 'SOA'}, 'Position', [465 9 100 20],...
       'FontSize', FontSize, 'Callback',@popup_system_Callback);
   htext_results = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'Results:',...
       'Position', [735 5 60 20], 'FontSize', FontSize);
   hpopup_results = uicontrol('Parent', hpanel_sim, 'Style','popupmenu',...
       'String', {'BER vs Transmitted Power', 'Power Penalty vs Modulator Cutoff Frequency',...
       'SOA/APD Gain'}, 'Position', [800 9 200 20],...
       'FontSize', FontSize, 'Callback', @popup_results_Callback);
   
   %% Block diagram  
   hpanel_blockdiagram = uipanel('Title', 'Block Diagram', 'FontSize', HeaderFontSize,...
       'Position', [0 1-hpanel_sim_height-.2 1 .2]);
   block_diagram = axes('Parent', hpanel_blockdiagram, 'Units', 'Pixels',...
       'Position', [5 5 fSize(1) 120], 'Color', 'w', 'Box', 'off');
   
   %% Results
   results = axes('Units', 'Pixels', 'Position', [fSize(1)-580 50 560 420], 'Box', 'on');
   
   %% Modulation Panel
   %% Transmitter Panel
   hpanel_tx = uipanel('Title', 'Transmitter', 'FontSize', HeaderFontSize,...
       'Position', [0 1-hpanel_sim_height-.2-.25 0.35 0.25]);
   
   [htext_lambda, hlambda] = table_entry(hpanel_tx, [0 0.9], 'Wavelength (nm):', 1310);
   [htext_fc, hfc] = table_entry(hpanel_tx, [0 0.75], 'Modulator Cutoff Frequency (GHz):', 30);
   [htext_kappa, hkappa] = table_entry(hpanel_tx, [0 0.6], 'Insertion Loss (dB):', 0);
   [hcheck_rex, htext_rex, hrex] = table_entry(hpanel_tx, [0 0.45], 'Extinction Ratio (dB):', -15, true);
   [hcheck_rin, htext_rin, hrin] = table_entry(hpanel_tx, [0 0.3], 'RIN (dB/Hz):', -150, true);
   [hcheck_chirp, htext_chirp, hchirp] = table_entry(hpanel_tx, [0 0.15], 'Modulator chirp (alpha):', -2);
   
   
  
   
   
% tx.alpha = 0; % chirp parameter
% tx.RIN = -150;  % dB/Hz
% tx.rexdB = -5;  % extinction ratio in dB. Defined as Pmin/Pmax
% 
% % Modulator frequency response
% tx.kappa = 1; % controls attenuation of I to P convertion
% tx.modulator.fc = 30e9; % modulator cut off frequency
   
%    align([hpanel_sim hpanel_blockdiagram], 'None', 'Distributed');
%    align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
   
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
 
    function varargout = table_entry(parent, norm_pos, str, value, checked)
        if nargin < 5
            checked = false;
        end
      checkbox = uicontrol('Parent', parent, 'Style', 'checkbox',...
       'String', '', 'Units', 'Normalized', 'Position', [norm_pos 0.05 0.1],...
       'FontSize', FontSize, 'Value', checked); % wavelength       
      text = uicontrol('Parent', parent, 'Style', 'text', 'HorizontalALignment', 'left',...
       'String', str, 'Units', 'Normalized', 'Position', [norm_pos(1)+0.07 norm_pos(2) 0.8 0.1],...
       'FontSize', FontSize); % wavelength
      edit = uicontrol('Parent', parent, 'Style', 'edit', 'String', num2str(value),...
           'HorizontalALignment', 'right', 'Units', 'Normalized', 'Position', [norm_pos(1)+0.85, norm_pos(2), 0.15 0.1],...
           'FontSize', FontSize); % wavelength
      align([checkbox, text, edit], 'None', 'Center')
        if nargout == 3
            varargout = {checkbox, text, edit};
        else
            varargout = {text, edit}
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