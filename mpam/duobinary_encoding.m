function z = duobinary_encoding(x, just_precoding)
%% duobinary encoding for multilevel signals
% x (PAM symbols) -> y (precoded duobinary symbols) -> z (duobinary
% filtering, if enabled)

% duobinary precoding
y = zeros(size(x));
for k = 2:length(x)
    y(k) = x(k)*sgn(y(k-1)) - y(k-1);
end

% duobinary filter
if exist('just_precoding', 'var') && just_precoding
    z = y;
else
    % duobinary filtering
    z = y + [0 y(1:end-1)];
end

% % duobinary precoding
% y = zeros(size(x));
% for k = 2:length(x)
%     y(k) = mod(x(k) - y(k-1), 4);
% %     y(k) = x(k) - mod(y(k-1), 4);
% end
% 
% % duobinary filter
% if exist('just_precoding', 'var') && just_precoding
%     z = y;
% else
%     % duobinary filtering
%     z = y + [0 y(1:end-1)];
% end
% 
% a = mod(z, 4);
% 1;

function s = sgn(x)
    % Differently from matlab's sign, this function
    % assigns a random -1, 1 if x == 0
    if x == 0
        s = [-1 1];
        s = s(randi([1 2], 1));
    else
        s = sign(x);
    end 
end
end