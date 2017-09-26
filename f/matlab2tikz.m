classdef matlab2tikz < handle
    properties
        xmin
        xmax
        ymin
        ymax
        xlabel
        ylabel
        legendLocation
        legendPosition = [];
        autoColorSelection = false;
        zplot = false
    end
    
    properties(GetAccess=private, Hidden)
        colorCount = 0;
    end        
    
    properties(Dependent, Hidden)
        tikzLegendLocation
    end
        
    properties(GetAccess=private, Constant)
        bl = 10; % break line character = \n
        tab = 9; % tab character = \t
        size = [4.52, 3.56]; % in  
        colorPalletNames = {'blue2', 'red2', 'green2', 'orange2',...
            'purple2', 'brown2', 'gray2', 'green3', 'black'}
        colorPalletValues = [51, 105, 232; 228,26,28; 0,153,37; 255,127,0;...
            152,78,163; 166,86,40; 153,153,153; 77,175,74; 0,0,0];  
    end
    
    properties (GetAccess=private)
        plots = [];
        annot = {};
    end
    
    methods
        function obj = matlab2tikz(h, zplot)
            %% Constructor
            if exist('zplot', 'var')
                obj.zplot = zplot;
            end
            
            if exist('h', 'var')
                obj = obj.extract(h);
            end 
        end
        
        function incColorCount(self)
            self.colorCount = mod(self.colorCount + 1, length(self.colorPalletNames));
        end
        
        function addheader(self, fileID)
            if self.zplot
                axisHeader = {'\begin{axis}[',...
                    'axis equal,',...
                    'axis lines*=middle,',...
                    'enlargelimits = false, clip=true,',...
                    sprintf('xmin=%.2f,', self.xmin),...
                    sprintf('xmax=%.2f,', self.xmax),...
                    sprintf('ymin=%.2f,', self.ymin),...
                    sprintf('ymax=%.2f,', self.ymax),...
                    'axis line style={->,>=stealth},',...
                    'xlabel={$\mathrm{Re}\{z\}$},',...
                    'ylabel={$\mathrm{Im}\{z\}$},',...
                    'every axis x label/.style={',...
                        'at={(ticklabel* cs:1)},',...
                        'anchor=north,',...
                    '},',...
                    'every axis y label/.style={',...
                        'at={(ticklabel* cs:1)},',...
                        'anchor=south,',...
                    '},',...
                    'xmajorgrids,','ymajorgrids,',...
                    'every outer y axis line/.append style={white!15!black},',...
                    'every y tick label/.append style={font=\color{white!15!black}},',...
                    'legend style={draw=white!15!black,fill=white,legend cell align=left}]'};
            else
                axisHeader = {'\begin{axis}[',...
                    sprintf('width=%.2fin,', self.size(1)),...
                    sprintf('height=%.2fin,', self.size(2)),...
                    'scale only axis,', 'separate axis lines,',...
                    'every outer x axis line/.append style={white!15!black},',...
                    'every x tick label/.append style={font=\color{white!15!black}},',...
                    sprintf('xmin=%.2f,', self.xmin),...
                    sprintf('xmax=%.2f,', self.xmax),...
                    sprintf('ymin=%.2f,', self.ymin),...
                    sprintf('ymax=%.2f,', self.ymax),...
                    ['xlabel={' self.xlabel '},'],...
                    ['ylabel={' self.ylabel '},'],...
                    'xmajorgrids,','ymajorgrids,',...
                    'every outer y axis line/.append style={white!15!black},',...
                    'every y tick label/.append style={font=\color{white!15!black}},',...                
                    'legend style={draw=white!15!black,fill=white,legend cell align=left}]'};

                    if ~isempty(self.legendPosition)
                        axisHeader{end} = ['legend style={draw=white!15!black,fill=white,legend cell align=left, at={(',...
                        sprintf('%f,%f', self.legendPosition(1), self.legendPosition(2)), ')},anchor=', 'south west', '}]']; % matlab seems to always anchor legend with south west
                    end
            end

            
            fprintf(fileID, '\\begin{tikzpicture}\n');
            for k = 1:length(axisHeader)
                fprintf(fileID, '%s\n', axisHeader{k});
            end
        end
              
        
        function self = setAxis(self, a)
            self.xmin = a(1);
            self.xmax = a(2);
            self.ymin = a(3);
            self.ymax = a(4);
        end
                
        function addplot(self, x, y, line, color, marker, label)
            %% Generate tikz table to be used in generating plots in latex
            % Inputs:
            % - x, y = data to be ploted
            % - line = matlab compatible type of line or symbol
            % i.e., {'-', '--', ':', 'none'}
            % - color = pgfplots compatible color {'blue', 'green', 'red', 'green2', 'orange',
            % 'purple', 'brown', 'gray'}. The colors passed will be matched
            % to colorPalletNames. First color of pallet is passed as
            % default if no match is found
            
            % - label = legend
            [color_rgb, color_name] = self.matchcolor(color);
            
            if exist('label', 'var')
                p = tikzplot(x, y, line, {color_rgb, color_name}, marker, label);          
            elseif exist('marker', 'var')
                p = tikzplot(x, y, line, {color_rgb, color_name}, marker);
            else
                p = tikzplot(x, y, line, {color_rgb, color_name});
            end
            
            self.plots = [self.plots; p];
        end
        
        function [color_rgb, color_name] = matchcolor(self, color)
            %% Match colors name to colors RGB
            if ~isempty(color) && ~strcmp(color, '')
                c = strfind(self.colorPalletNames, color);
                for k = 1:length(c)
                    if ~isempty(c{k})
                        color_rgb = self.colorPalletValues(k, :);
                        break
                    end
                end
                color_name = self.colorPalletNames{k};
            else
                color_rgb = self.colorPalletValues(self.colorCount+1, :);
                color_name = self.colorPalletNames{self.colorCount+1};
                self.incColorCount();
            end
            color_rgb = color_rgb/255;
        end
        
        function self = extract(self, h, type)
            %% extract properties from handle h
            if ~exist('type', 'var')
                type = 'curves';
            end
            
            if ~strcmp(type, 'just axis')
                children = allchild(h);
                for k = 1:length(children)
                    if strcmp(get(children(k), 'type'), 'line')
                        if self.zplot && strcmpi(get(children(k), 'LineStyle'), ':')
                            continue % ignore unit circle
                        end
                        p = tikzplot(children(k));
                        if self.zplot && strcmpi(get(children(k), 'Marker'), 'x')
                            p.label = '__poles__';
                        elseif self.zplot && strcmpi(get(children(k), 'Marker'), 'o')
                            p.label = '__zeros__';
                        end
                        
                        if self.autoColorSelection % select colors automatically
                            p.color = self.colorPalletValues(self.colorCount+1, :)/255;
                            p.colorName = self.colorPalletNames{self.colorCount+1};
                            self.colorCount = mod(self.colorCount + 1, length(self.colorPalletNames));
                        else % uses matlab's colors
                            p.color = get(children(k), 'Color');
                            p.colorName = ['matlabColor' num2str(k)];                            
                        end
                        self.plots = [self.plots; p];                    
                    elseif strcmpi(get(children(k), 'type'), 'text')
                        a = self.extract_annot(children(k));
                        self.annot = [self.annot, a];
                    end
                end  
            end
                       
            % xmax, xmin
            x = get(h, 'XLim');
            self.xmin = x(1);
            self.xmax = x(2);
            % ymax, ymin
            y = get(h, 'YLim');
            self.ymin = y(1);
            self.ymax = y(2);
            % 
            self.xlabel = get(get(h, 'xlabel'), 'String');
            self.ylabel = get(get(h, 'ylabel'), 'String');
            % 
            self.legendLocation = get(legend(h), 'Location');
            self.legendPosition = get(legend(h), 'Position');
        end
        
        function a = extract_annot(~, h)
            a = sprintf('\\node[align=%s, anchor=south] at(axis cs: %s, %s) {\\scriptsize $%s$};',...
                h.HorizontalAlignment,...
                num2str(h.Position(1)),...
                num2str(h.Position(2)),...
                h.String);            
        end
        
        function clearPlots(self)
            self.plots = [];
            self.colorCount = 0;
        end            
        
        function writeplots(self, fileID)
            if self.zplot
                fprintf(fileID, '\\draw (axis cs:0,0) circle [black!50, dashed, line width=2pt, radius=1];\n');
            end
            
            for n = 1:length(self.plots)
                p = self.plots(n);
                fprintf(fileID, '%s\n', p.addplot());
            end
        end
        
        function writeannot(self, fileID)
            if isempty(self.annot)
                return;
            end
            
            fprintf(fileID, '%% Annotations\n');
            for n = 1:length(self.annot)
                fprintf(fileID, '%s\n', self.annot{n});
            end
        end
        
        function loc = get.tikzLegendLocation(self)
            m = {'northwest', 'southwest', 'southeast', 'northeast'};
            t = {'north west', 'south west', 'south east', 'north east'};
            loc = t{strcmp(self.legendLocation, m)};
            if isempty(loc)
                loc = 'north west';
            end
        end
        
        function write(self, filename)            
            try
                fileID = fopen(filename, 'w');
            catch e
                warning('matlab2tikz/write: Error opening the file %f:\n%s\n', filename, e.message)
                return
            end

            self.addheader(fileID)
            self.writeplots(fileID);
            self.writeannot(fileID);
            
            fprintf(fileID, '\\end{axis}\n\\end{tikzpicture}');
            fclose(fileID);
        end
        
        function write_tables(self, filename, param) % filename without extension
            if exist('param', 'var') && strcmpi(param, 'same x')
                try
                   filename_complete = [filename '.dat'];
                   fileID = fopen(filename_complete, 'w');
               catch e
                    warning('matlab2tikz/write_table: Error opening the file %f:\n%s\n', filename, e.message)
                    fclose(fileID)
                    return
                end
                
                x = self.plots(1).x;
                str = 'x';
                for n = 1:length(self.plots)
                    str = [str self.tab self.plots(n).label];
                end
                str = [str self.bl];
                
                for k = 1:length(x)
                    str = [str num2str(x(k))];
                    for n = 1:length(self.plots)
                        str = [str self.tab num2str(self.plots(n).y(k))];
                    end
                    str = [str self.bl];
                end 
                
                fprintf(fileID, '%s', str);
                fclose(fileID);
                return
            end
            
            for n = 1:length(self.plots)
                try
                    if strcmpi(self.plots(n).label, '')
                        filename_complete = [filename '.dat'];
                    else
                        lab = strrep(self.plots(n).label, ' ',''); % remove spaces
                        lab = strrep(lab, '-',''); % remove hyphens
                        lab = strrep(lab, '.',''); % remove periods
                        filename_complete = [filename '_' lower(lab) '.dat'];
                    end
                    fileID = fopen(filename_complete, 'w');

                    p = self.plots(n);

                    % Add data table
                    str = '';
                    for k = 1:length(p.x)
                        str = [str self.tab num2str(p.x(k)) self.tab num2str(p.y(k)) self.bl];
                    end
                    fprintf(fileID, '%s', str);
                    fclose(fileID);
                catch e
                    warning('matlab2tikz/write_table: Error opening the file %f:\n%s\n', filename, e.message)
                    fclose(fileID)
                    return
                end
            end
        end
    end
end