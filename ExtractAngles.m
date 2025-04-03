if not(isfolder("NoPain_preprocessed"))
mkdir("NoPain_preprocessed")
end

if not(isfolder("Pain_preprocessed"))
mkdir("Pain_preprocessed")
end

readCSV('NoPain');
readCSV('Pain');

function readCSV(patientsgroup)

rootfolder = pwd;
patients= dir(fullfile(rootfolder,patientsgroup,'*.csv'));

for i = 1:length(patients)
    path = fullfile(rootfolder,patientsgroup);
    opts = detectImportOptions(fullfile(path, patients(i).name));
    opts = setvartype(opts,'string');
    patient_data = readtable(fullfile(path, patients(i).name),opts);
    
    %Search for the first apple that appears in each patient
    for j= 1:size(patient_data)
        if contains(patient_data(j,:).Var3, "spawn,Apple 0")
            patient_data_initiated= patient_data(j:height(patient_data),:);
        end
    end

    %Search for the headset (H) data and location of stored apples
    search = ["H","stored"];
    headset_data = patient_data_initiated(contains(patient_data_initiated.Var3,search), :);
    % Extract the Euler angles (x, y, z) from the headset data and the positition of the apples
    euler_angles = headset_data(:, [1 2 3 11:13]);
    % We save the extracted data to a new CSV file
    cd (rootfolder+"\"+patientsgroup+"_preprocessed\");
    filename= strcat(erase(patients(i).name,".csv"),'_preprocessed.csv');
    writetable(euler_angles, filename);
end
cd (rootfolder);

end
