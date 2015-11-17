classdef matlab2tikz < handle
    properties
        xmin
        xmax
        ymin
        ymax
        xlabel
        ylabel
        legendLocation
    end
    
    properties(GetAccess=private, Constant)
        size = [4.52, 3.56]; % in  
        colorPalletNames = {'blue', 'green', 'red', 'green2', 'orange', 'purple',...
            'brown', 'gray'}
        colorPalletValues = [51, 105, 232; 0, 153, 37; 228,26,28; 77,175,74;...
            255,127,0; 152,78,163; 166,86,40; 153,153,153];  
    end
    properties (GetAccess=private)
        plots = {}
    end
    
    methods
        function obj = matlab2tikz()
            %% Constructor
            
        end
        
        function addheader(self, fileID)
            axisHeader = {'\begin{axis}[',...
                sprint('width=%.2fin', self.size(1)),...
                sprint('height=%.2fin', self.size(2)),...
                'scale only axis', 'separate axis lines',...
                'every outer x axis line/.append style={white!15!black}',...
                'every x tick label/.append style={font=\color{white!15!black}}',...
                sprint('xmin=%.2fin', self.xmin),...
                sprint('xmax=%.2fin', self.xmax),...
                sprint('ymin=%.2fin', self.ymin),...
                sprint('ymax=%.2fin', self.ymax),...
                ['xlabel={' self.xlabel '}'],...
                ['ylabel={' self.ylabel '}'],...
                'xmajorgrids','ymajorgrids',...
                'every outer y axis line/.append style={white!15!black}',...
                'every y tick label/.append style={font=\color{white!15!black}}',...                
                'legend style={draw=white!15!black,fill=white,legend cell align=left, at={(0.03,0.97)},anchor=north west}',...
                ']'};
            
            fprintf(fileID, '\begin{tikzpicture}\n');
            for k = 1:length(axisHeader)
                fprintf(fileID, '\t%s\n', axisHeader{k});
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
            % - line = pgfplots compatible type of line or symbol
            % i.e., {'solid', 'dashed', 'dotted', 'none', circle, square*, etc}
            % - color = {'blue', 'green', 'red', 'green2', 'orange',
            % 'purple', 'brown', 'gray'} or any other native pgfplots color
            % - label = legend
            
            self.plots = [self.plots; tikzplot(x, y, line, color, marker, label)];
        end
        
        function writeplots(self, fileID)
            for n = 1:length(self.plots)
                p = self.plots{n};
                fprintf(fileID, '\t%s\n', p.addplot());
            end
        end
        
        function write(self, filename)
            try
                fileID = fopen(filename, 'w');
            catch e
                warning('matlab2tikz/write: Error opening the file %f:\n%s\n', filename, e.message)
                return
            end

            fprintf(fileID, 'table[row sep=crcr]{\n');

            self.addheader(fileID)
            self.writeplots(fileID);
            
            fprintf(fileID, '};\n\end{tikzpicture}');
            fclose(fileID);
        end
        
    end
end