function sim_100G_single_laser2
   %  Create and then hide the GUI as it is being constructed.
   f = figure('Visible','off','Position',[360,500,900,450]);
 
   HeaderFontSize = 12;
   FontSize = 11;
   
   % Simulation
   hpanel_sim_height = .11;
   hpanel_sim = uipanel('Title', 'Simulation', 'FontSize', HeaderFontSize, 'Position',[0 1-hpanel_sim_height 1 hpanel_sim_height]);
   htext_modulation = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'Modulation Type:', 'Position', [0 5 120 20], 'FontSize', FontSize);
   hpopup_modulation = uicontrol('Parent', hpanel_sim, 'Style','popupmenu',...
       'String', {'M-PAM', 'M-CAP', 'DMT/OFDM'}, 'Position', [125 9 100 20],...
       'FontSize', FontSize, 'Callback',@popup_modulation_Callback);
   htext_system = uicontrol('Parent', hpanel_sim, 'Style', 'text', 'String', 'System:', 'Position', [150 5 60 20], 'FontSize', FontSize);
   hpopup_system = uicontrol('Parent', hpanel_sim, 'Style','popupmenu',...
       'String', {'Basic', 'APD', 'SOA'}, 'Position', [210 9 100 20],...
       'FontSize', FontSize, 'Callback',@popup_modulation_Callback);
   
   %  Construct the components.
   hsurf = uicontrol('Style','pushbutton','String','Surf',...
          'Position',[315,220,70,25],...
          'Callback',@surfbutton_Callback);
   hmesh = uicontrol('Style','pushbutton','String','Mesh',...
          'Position',[315,180,70,25],...
          'Callback',@meshbutton_Callback);
   hcontour = uicontrol('Style','pushbutton',...
          'String','Countour',...
          'Position',[315,135,70,25],...
          'Callback',@contourbutton_Callback); 
   htext = uicontrol('Style','text','String','Select Data',...
          'Position',[325,90,60,15]);
   hpopup = uicontrol('Style','popupmenu',...
          'String',{'Peaks','Membrane','Sinc'},...
          'Position',[300,50,100,25],...
          'Callback',@popup_menu_Callback);
   ha = axes('Units','Pixels','Position',[50,60,200,185]); 
   align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
   
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
   current_data = peaks_data;
   surf(current_data);
   % Assign the GUI a name to appear in the window title.
   f.Name = 'Simple GUI';
   % Move the GUI to the center of the screen.
   movegui(f,'center')
   % Make the GUI visible.
   f.Visible = 'on';
 
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