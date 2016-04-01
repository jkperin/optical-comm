%% Positive-negative photodiode
classdef pin < apd % inherits methods and properties from apd
    %% Positive-negative photodiode inherits methods and properties from class apd
    properties
        % defined in apd.m
        % R    % responsivity
        % Id   % dark current  
        % BW   % bandwidth
    end
    methods
        function self = pin(R, Id, BW)
            if ~exist('R', 'var')
                R = 1;
            end
            
            if ~exist('BW', 'var')
                BW = Inf;
            end
            
            if ~exist('Id', 'var')
                Id = 0;
            end
            % Call constructor from APD
            % apd(GaindB, ka, BW, R, Id)
            self@apd(0, 0, BW, R, Id);
        end       
    end
end