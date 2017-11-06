function improve_plot(h, other)

children = allchild(h);
for k = 1:length(children)
    children(k).LineWidth = 2;
end

set(h, 'box', 'on')
set(h, 'FontSize', 12)
xl = get(h, 'xlabel');
xl.FontSize = 12;
yl = get(h, 'ylabel');
yl.FontSize = 12;

if exist('other', 'var')
    if iscell(other)
        for k = 1:length(other)
            eval_cmd(other{k})
        end
    else
         eval_cmd(other)
    end
end

function eval_cmd(cmd)
    try
        eval(cmd)
    catch e
        warning('Error %s while evaluating command %s\n', e.message, cmd)
    end
end
end
    
