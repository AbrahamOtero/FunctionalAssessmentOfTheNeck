global B_1_start B_1_end B_2_start B_2_end B_3_start B_3_end 

B_1_start=0;
B_1_end = 0.1;
B_2_start = 0.1;
B_2_end = 0.5;
B_3_start = 0.5;
B_3_end = 5;

average_NoPain= calculateSpectralPower('NoPain');
average_Pain= calculateSpectralPower('Pain');
[hypothesis_values_freq, p_values_freq,averagePower,stdPower]= test_frequency(average_NoPain, average_Pain);

duration_NoPain= calculateDuration('NoPain');
duration_Pain= calculateDuration('Pain');
[hypothesis_value_duration, p_value_duration]= test_duration(duration_NoPain,duration_Pain);

%Correct p-values
p_values_freq_adj = pval_adjust(p_values_freq, 'fdr');

disp('P values (Frequecy)')
p_values_freq
disp('Corrected P values (Frquency)')
p_values_freq_adj

disp('P value (Time)')
p_value_duration
            

function avgpower= calculateSpectralPower (groupname)

global B_1_start B_1_end B_2_start B_2_end B_3_start B_3_end

window_length =512;
overlap = 128;
f = [B_1_start:0.01:B_1_end, B_2_start:0.01:B_2_end, B_3_start:0.01:B_3_end];

avgpower=zeros(30,3,4);
patient_index=1;

rootfolder = pwd;
addpath(rootfolder);
groupfolder= rootfolder+"\"+groupname+"_preprocessed\";
cd (groupfolder);

for i=1:60
    if isfile("Pte_"+i+"_preprocessed.csv")
        read_data= readtable("Pte_"+i+"_preprocessed.csv");
        read_data.Var3 = string(read_data.Var3);
        read_data.Var1 = num2str(read_data.Var1);
        read_data_array = table2array(read_data);
        rowsEuler= ~contains(read_data_array(:, 3),"stored");
        patient_xyz= str2double(strrep(read_data_array(rowsEuler,:),',','.')); 

        %We discard the first 30 seconds to ignore quiet at the beginning
        fs= 50; 
        initial_30_seconds = 30 * fs;
        patient_x=patient_xyz(:,4)-patient_xyz(1,4);
        patient_x = patient_x(initial_30_seconds+1:end);
        patient_y=patient_xyz(:,5)-patient_xyz(1,5);
        patient_y = patient_y(initial_30_seconds+1:end);
        patient_z=patient_xyz(:,6)-patient_xyz(1,6);
        patient_z = patient_z(initial_30_seconds+1:end);
        %Avoid +180 -180 oscilations
        patient_x  = wrapTo180(patient_x);
        patient_y  = wrapTo180(patient_y);
        patient_z  = wrapTo180(patient_z);
        %Total movement
        patient_xyz= sqrt(patient_x.^2 + patient_y.^2 + patient_z.^2);
        %Time series are centred
        patient_x_centred=patient_x-mean(patient_x);
        patient_y_centred=patient_y-mean(patient_y);
        patient_z_centred=patient_z-mean(patient_z);
        patient_xyz_centred=patient_xyz-mean(patient_xyz);

        [S, F, T] = spectrogram(patient_x_centred,window_length,overlap,f,fs);
        avgpower(patient_index,:,1) = [mean(sum(abs(S(F >= B_1_start & F <= B_1_end, :)), 1)), mean(sum(abs(S(F >B_2_start & F <= B_2_end, :)), 1)), mean(sum(abs(S(F > B_3_start & F <= B_3_end, :)), 1))];

        [S, F, T] = spectrogram(patient_y_centred,window_length,overlap,f,fs);
        avgpower(patient_index,:,2) = [mean(sum(abs(S(F >= B_1_start & F <= B_1_end, :)), 1)), mean(sum(abs(S(F >B_2_start & F <= B_2_end, :)), 1)), mean(sum(abs(S(F > B_3_start & F <= B_3_end, :)), 1))];

        [S, F, T] = spectrogram(patient_z_centred,window_length,overlap,f,fs);
        avgpower(patient_index,:,3) = [mean(sum(abs(S(F >= B_1_start & F <= B_1_end, :)), 1)), mean(sum(abs(S(F >B_2_start & F <= B_2_end, :)), 1)), mean(sum(abs(S(F > B_3_start & F <= B_3_end, :)), 1))];

        [S, F, T] = spectrogram(patient_xyz_centred,window_length,overlap,f,fs);
        avgpower(patient_index,:,4) = [mean(sum(abs(S(F >= B_1_start & F <= B_1_end, :)), 1)), mean(sum(abs(S(F >B_2_start & F <= B_2_end, :)), 1)), mean(sum(abs(S(F > B_3_start & F <= B_3_end, :)), 1))];

        patient_index=patient_index+1;
    end 
end
cd (rootfolder)
end 


function [hypothesis_values_freq,p_values_freq,averagePower,stdPower]= test_frequency(avg_NoPain, avg_Pain)

h_NoPain= zeros(4,3);
h_Pain= zeros(4,3);
hypothesis_values_freq= zeros(4,3);
p_values_freq= zeros(4,3);
averagePower=zeros(2,3,4);
stdPower=zeros(2,3,4);

%Loop through each axis
for k = 1:4 
    %Loop through each frequency band
    for j = 1:3 
        %Check normality
        h_NoPain(k,j) = adtest(avg_NoPain(:,j,k));
        h_Pain(k,j) = adtest(avg_Pain(:,j,k));

        %T-test or Whitney-U-test based on normality
        if h_NoPain(k,j)|| h_Pain(k,j)
            [p_values_freq(k,j),hypothesis_values_freq(k,j)]= ranksum(avg_NoPain(:,j,k), avg_Pain(:,j,k));
            %fprintf('k: %i, j: %.i not normal', k, j);
        else 
            [hypothesis_values_freq(k,j),p_values_freq(k,j)]= ttest2(avg_NoPain(:,j,k), avg_Pain(:,j,k));
        end 
        
        %We calculate mean and std deviation 
        averagePower(1,j,k)= mean(avg_NoPain(:,j,k));
        averagePower(2,j,k)= mean(avg_Pain(:,j,k));
        stdPower(1,j,k)= std(avg_NoPain(:,j,k));
        stdPower(2,j,k)= std(avg_Pain(:,j,k));       
    end
end
end


%Duration of the spectrograms 
function duration_patients= calculateDuration(groupname)

rootfolder = pwd;
addpath(rootfolder);
groupfolder= rootfolder+"\"+groupname+"_preprocessed\";
cd (groupfolder);

duration_patients= zeros(1,30);
index=1;

for i=1:60
    if isfile("Pte_"+i+"_preprocessed.csv")
        read_data= readtable("Pte_"+i+"_preprocessed.csv");
        read_data.Var3 = string(read_data.Var3);
        read_data.Var1 = num2str(read_data.Var1);
        read_data_array = table2array(read_data);
        patient_xyz= str2double(strrep(read_data_array,',','.')); 
       
        duration_patients(index)= (patient_xyz(size(patient_xyz,1),2)-patient_xyz(1,2))/1000;
        index= index+1;
    end 
end
cd (rootfolder)
end 


function [hypothesis_value_duration, p_value_duration]= test_duration(duration_NoPain,duration_Pain)
%Check normality
h_NoPain = adtest(duration_NoPain);
h_Pain = adtest(duration_Pain);

%Perform t-test or Whitney-U-test based on normality assumptions
if h_NoPain|| h_Pain
    [p_value_duration,hypothesis_value_duration]= ranksum(duration_NoPain,duration_Pain);
    %disp('not noraml')
else 
    [hypothesis_value_duration,p_value_duration]= ttest2(duration_NoPain,duration_Pain);
end      
end
   