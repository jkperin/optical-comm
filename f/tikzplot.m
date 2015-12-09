classdef tikzplot
    properties
        x
        y
        color = [0 0 0]; % rgb values
        colorName = 'black'; % must be tikz compatible (not 'b', 'r')
        marker = 'none' % must be tikz compatible (not 's', 'v')
        markerSize = 2 % in pt
        line = 'solid' % must be tikz compatible (not '-')
        lineWidth = 1.5 % in pt
        label = '';
    end
    
    properties (Constant, Hidden)
        bl = 10;
        tab = 9;
        mlineStyles = {'none', '-', '--', ':'};
        tikzlineStyles = {'none', 'solid', 'dashed', 'dotted'};
        mmarkerStyles = {'none', 'o', 's'};
        tikzmarkerStyles = {'none', '*', 'square*'}; % filled by default
    end
    
    methods
        function obj = tikzplot(x, y, line, color, marker, label)
            %% Constructor
            if nargin == 1
                obj = obj.extract(x);
                return
            end
            
            obj.x = x;
            obj.y = y;
            assert(~(isempty(x) || isempty(y) || length(x) ~= length(y)),...
                'x and y dimensions must be non-zero and must have the same length')
            
            if exist('line', 'var')
                obj.line = obj.m2tikzline(line);
            end
            if exist('color', 'var')
                obj.color = color{1};
                obj.colorName = color{2};
            end
            if exist('marker', 'var')
                obj.marker = obj.m2tikzmark(marker);
            end
            if exist('label', 'var')
                obj.label = label;           
            end
        end
        
        function tikzline = m2tikzline(self, mline)
            ind = strcmp(self.mlineStyles, mline);
            if isempty(ind)
                ind = 2; % assume solid by default
            end
            tikzline = self.tikzlineStyles{ind};
        end
        
        function tikzmark = m2tikzmark(self, mmark)
            ind = strcmp(self.mmarkerStyles, mmark);
            if isempty(ind)
                ind = 1; % assume none by default
            end
            tikzmark = self.tikzmarkerStyles{ind};
        end
                
        function self = extract(self, h)
            % LineWidth is not imported 
            self.color = get(h, 'Color');
            self.line = self.m2tikzline(get(h, 'LineStyle'));
            self.marker =  self.m2tikzmark(get(h, 'Marker')); 
%             self.markerSize = get(h, 'MarkerSize');
            self.x = get(h, 'XData');
            self.y = get(h, 'Ydata');
            self.label = get(h, 'DisplayName');
            
            if isempty(self.line)
                self.line = 'solid';
            end
            if isempty(self.marker)
                self.marker = 'none';
            end
        end
        
        function str = addplot(self)
            %% Generate tikz table to be used in generating plots in latex
          
            % Color
%             str = ['\addplot [color=' self.color];
            str = sprintf('\\definecolor{%s}{rgb}{%f,%f,%f}\n',self.colorName,...
                self.color(1), self.color(2), self.color(3));
            str = [str '\addplot [color=' self.colorName];
            
            % Line type
            if strcmp(self.line, 'none')
                str = [str ', only marks'];
            else
                str = [str ', ', self.line];
            end
            
            % Line width
            str = [str ', line width=', num2str(self.lineWidth), 'pt'];
            
            % Marker
            if ~strcmp(self.marker, 'none')
                str = [str ', mark=', self.marker, ', mark size=', num2str(self.markerSize), 'pt'];
            end                  
            
            % maker marker line solid
            if strcmp(self.line, 'dashed')
                str = [str ', mark options={solid}'];
            end
            
            % If no legend, then forget plot
            if strcmp(self.label, '')
                str = [str ', forget plot'];
            end            
            
            % Close bracket
            str = [str ']' self.bl];
            
            % Add data table
            str = [str, 'table[row sep=crcr]{', self.bl];
            for k = 1:length(self.x)
                str = [str self.tab num2str(self.x(k)) ' ' num2str(self.y(k)) ' \\', self.bl];
            end
            
            % Add legend
            if ~strcmp(self.label, '')
                str = [str '}; \addlegendentry{' self.label '}', self.bl];
            else
                str = [str '};', self.bl];
            end
        end
    end
end