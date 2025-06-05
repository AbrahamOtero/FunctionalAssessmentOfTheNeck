
function myPlotSpectrogram (s,f,t, label)
    fontSize = 16;
    surf(t, f, 20*log10(abs(s)), 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)', 'FontSize', fontSize);
    colorbar;
    ylabel(label, 'FontSize', fontSize);
    clim([-1, 80])
    fontSize = 13; 
    set(gca, 'FontSize', fontSize);
end