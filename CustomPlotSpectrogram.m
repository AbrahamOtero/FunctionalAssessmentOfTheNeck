
function myPlotSpectrogram (s,f,t, label)
    surf(t, f, 20*log10(abs(s)), 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    colorbar;
    ylabel(label);
    clim([-10, 90])

end