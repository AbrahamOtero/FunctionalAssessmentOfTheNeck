plotSpectrograms('NoPain');
plotSpectrograms('Pain');

function plotSpectrograms (groupname)
rootfolder = pwd;
addpath(rootfolder);
groupfolder= rootfolder+"\"+groupname+"_preprocessed\";
cd (groupfolder);

if not(isfolder(groupname+"_spectrograms"))
mkdir(groupname+"_spectrograms")
end

global B_1_start B_1_end B_2_start B_2_end B_3_start B_3_end

f = [B_1_start:0.01:B_1_end, B_2_start:0.01:B_2_end, B_3_start:0.01:B_3_end];

for i=1:60
    if isfile("Pte_"+i+"_preprocessed.csv")
        read_data= readtable("Pte_"+i+"_preprocessed.csv");

        read_data.Var3 = string(read_data.Var3);
        read_data.Var1 = num2str(read_data.Var1);
        read_data_array = table2array(read_data);
        rowsEuler= ~contains(read_data_array(:, 3), "stored");
        rowsApples= contains(read_data_array(:, 3), "stored");
        patient_xyz= str2double(strrep(read_data_array(rowsEuler,:),',','.')); 
        samples_storedApples= str2double(strrep(read_data_array(rowsApples,1:2),',','.'));

        fs= 50;
        window_length =128;
        overlap = 18;
        figure('Position', [100, 100, 800, 600]);
        %figure;
        sub1=subplot(4,1,1);
        patient_x=patient_xyz(:,4)-patient_xyz(1,4);
        [s,f,t] =spectrogram(patient_x,window_length,overlap,f,  fs, 'yaxis');
        CustomPlotSpectrogram(s,f,t, 'X axis (Hz)');
        time_storedApples =(samples_storedApples(:, 2)-(patient_xyz(1,2)))/1000; 
        y_lim= get(gca, 'YLim');
        z_lim= max(max(20*log10(abs(s))));
        sub1.Tag='1';
        ax1 = findobj(gcf,'Tag','1'); 
        hold(ax1,'on');
        for index = 1:length(time_storedApples)
        [~, closest_index] = min(abs(t - time_storedApples(index)));
        line([t(closest_index),t(closest_index)], y_lim, [z_lim z_lim],'Color', 'w', 'LineWidth', 1.5);
        end

        sub2= subplot(4,1,2);
        patient_y=patient_xyz(:,5)-patient_xyz(1,5);
        [s,f,t] = spectrogram(patient_y,window_length,overlap, f, fs,'yaxis');
        CustomPlotSpectrogram(s,f,t, 'Y axis (Hz)');      
        sub2.Tag='1';
        ax2 = findobj(gcf,'Tag','1'); 
        hold(ax2,'on');
        for index = 1:length(time_storedApples)
        [~, closest_index] = min(abs(t - time_storedApples(index)));
        line([t(closest_index),t(closest_index)], y_lim, [z_lim z_lim],'Color', 'w', 'LineWidth', 1.5);
        end
        
        
        sub3= subplot(4,1,3);
        patient_z=patient_xyz(:,6)-patient_xyz(1,6);
        [s,f,t] = spectrogram(patient_z,window_length,overlap, f, fs,'yaxis');
        CustomPlotSpectrogram(s,f,t, 'Z axis(Hz)');
        sub3.Tag='1';
        ax3 = findobj(gcf,'Tag','1'); 
        hold(ax3,'on');
        for index = 1:length(time_storedApples)
        [~, closest_index] = min(abs(t - time_storedApples(index)));
        line([t(closest_index),t(closest_index)], y_lim, [z_lim z_lim],'Color', 'w', 'LineWidth', 1.5);
        end

        sub4= subplot(4,1,4);
        [s,f,t] = spectrogram((abs(patient_x)+abs(patient_y)+abs(patient_z)),window_length,overlap, f, fs,'yaxis');
        CustomPlotSpectrogram(s,f,t, 'Total (Hz)');
        sub4.Tag='1';
        ax4 = findobj(gcf,'Tag','1'); 
        hold(ax4,'on');
         for index = 1:length(time_storedApples)
        [~, closest_index] = min(abs(t - time_storedApples(index)));
        line([t(closest_index),t(closest_index)], y_lim, [z_lim z_lim],'Color', 'w', 'LineWidth', 1.5);
        end

        cd (groupfolder+"\"+groupname+"_spectrograms\");
        saveas(gcf,"Spectrogram_Pte"+i+".png");
        close(gcf);
        cd (groupfolder);
    end 
end
cd (rootfolder)

end 
