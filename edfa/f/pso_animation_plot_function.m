function stop = pso_animation_plot_function(optimValues, state, lnm)

    stop = false;
    iteration = optimValues.iteration;
    swarm = optimValues.swarm;
    swarmfvals = abs(optimValues.swarmfvals);
    
    if mod(iteration, 50) ~= 0
        return;
    end
        
    
    figure(1), set(gcf, 'color', 'w'), hold off, box on
    plot(lnm, swarm)
    xlabel('Wavelength (nm)', 'FontSize', 14)
    ylabel('Input signal power (dBm)', 'FontSize', 14)
    set(gca, 'FontSize', 14)
    axis([lnm(1) lnm(end) -30 -10])
    drawnow
    
    
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    outfile = 'pso.gif';
     
    % On the first loop, create the file. In subsequent loops, append.
    if exist('pso.gif', 'file') == 0
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end
end
