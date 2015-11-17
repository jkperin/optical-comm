classdef tikzplot
    properties
        x
        y
        color = 'black'
        marker = 'none'
        markerSize = 2 % in pt
        line = 'solid'
        lineWidth = 1.5 % in pt
        label = '';
    end
    
    methods
        function obj = tikzplot(x, y, line, color, marker, label)
            %% Constructor
            obj.x = x;
            obj.y = y;
            assert(isempty(x) || isempty(y) || length(x) ~= length(y), 'x and y dimensions must be non-zero and must have the same length')
            
            if exist('line', 'var')
                obj.line = line;
            end
            if exist('color', 'var')
                obj.color = color;
            end
            if exit('marker', 'var')
                obj.marker = marker;
            end
            if exist('label', 'var')
                obj.label = label;           
            end
        end
        
        function str = addplot(self)
            %% Generate tikz table to be used in generating plots in latex
            % Inputs:
            % - x, y = data to be ploted
            % - line = pgfplots compatible type of line or symbol
            % i.e., {'solid', 'dashed', 'dotted', 'none', circle, square*, etc}
            % - color = {'blue', 'green', 'red', 'green2', 'orange',
            % 'purple', 'brown', 'gray'} or any other native pgfplots color
            % - label = legend
            
            str = ['\addplot [color=' self.color];
            % Line type
            if self.line == 'none'
                str = [str ', only marks'];
            else
                str = [str ', ', self.line];
            end
            
            % Line width
            str = [str ', line width=', num2str(self.lineWidth), 'pt'];
            
            % Marker
            if self.mark ~= 'none'
                str = [str ', mark=', self.mark, ', mark size=', num2str(self.markerSize), 'pt'];
            end                      
            str = [str ']\n'];
            
            % Add data table
            str = [str, 'table[row sep=crcr]{\n'];
            for k = 1:length(self.x)
                str = [str '\t' num2str(self.x(k)) ' ' num2str(self.y(k)) ' \\\\\n'];
            end
            
            % Add legend
            if self.label == '' 
                str = [str '}; \addlegendentry{' self.label '}\n'];
            else
                str = [str '};\n'];
            end
        end
    end
end