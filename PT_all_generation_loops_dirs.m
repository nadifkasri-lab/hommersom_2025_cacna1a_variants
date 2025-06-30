%This code generates .csv files with spike timings from .spk files as generated for Axion
%Biosystems MEA reocrdings. It is based on AxionFileLoader, which is available through
%https://github.com/axionbio/AxionFileLoader

top_level_dir = 'insert/folder/name'; % This is the higher-level parent folder
cd(top_level_dir)

plate_dirs = dir(top_level_dir); 
plate_dirs = plate_dirs([plate_dirs.isdir]);
plate_dirs = plate_dirs(~ismember({plate_dirs.name}, {'.', '..'})); % Remove '.' and '..'

% Loop over each "plate" folder
for p = 1:length(plate_dirs)
    plate_num = plate_dirs(p).name; % Now the folder name is the plate number
    D = fullfile(top_level_dir, plate_num);
    cd(D)
    
    files = dir(D); 
    dirFlags = [files.isdir];
    subDirs = files(dirFlags);
    folders = subDirs(~ismember({subDirs.name}, {'.','..'})); % Actual subfolders
    L = length(folders);
    
    for f = 1:L
        spk_dir = folders(f).name;
        spk_path = fullfile(D, spk_dir, plate_num);
        
        % Compatibility-safe directory check
        if exist(spk_path, 'dir') ~= 7
            fprintf('Skipping missing directory: %s\n', spk_path);
            continue
        end
        
        cd(spk_path)
        spk_files = dir('*.spk');
        S = length(spk_files);
        
        for spk = 1:S
            spk_name = spk_files(spk).name;
            SpkPath = fullfile(spk_path, spk_name);
            
            AllData = AxisFile(SpkPath).SpikeData.LoadData;
            final_results = [];
            [nwr, nwc, nec, ner] = size(AllData);
            
            for wr = 1:nwr
                for wc = 1:nwc
                    for ec = 1:nec
                        for er = 1:ner
                            data = AllData{wr, wc, ec, er};
                            if ~isempty(data)
                                [t, v] = data.GetTimeVoltageVector;
                                timestamp = t(1, :);
                                timestamp_length = length(timestamp);
                                Channel_Label = str2double(strcat(num2str(ec), num2str(er)));
                                Channel_Label = repelem(Channel_Label, timestamp_length);
                                Well_Label = str2double(strcat(num2str(wr), num2str(wc)));
                                Well_Label = repelem(Well_Label, timestamp_length);
                                min_amplitude = min(v);
                                max_amplitude = max(v);
                                peak_to_peak_amplitude = max_amplitude - min_amplitude;
                                combined_data = transpose([Channel_Label; Well_Label; timestamp; max_amplitude; min_amplitude; peak_to_peak_amplitude]);
                                final_results = [final_results; combined_data];
                            else
                                fprintf('no spikes\n');
                            end
                        end
                    end
                end
            end
            
            final_results = array2table(final_results);
            final_results.Properties.VariableNames = {'Channel_Label','Well_Label','Timestamp','Maximum_Amplitude','Minimum_Amplitude','Peak_to_peak_Amplitude'};
            
            % Make subdirectory and write result
            ix = strfind(spk_name, '_');
            if length(ix) >= 2
                l = length(spk_name);
                dir_name = spk_name((ix(end-1)+1):(l-4));
            else
                dir_name = 'Default';
            end
            output_dir = fullfile(spk_path, dir_name);
            if exist(output_dir, 'dir') ~= 7
                mkdir(output_dir)
            end
            table_name = fullfile(output_dir, 'Axion_PT_all.csv');
            writetable(final_results, table_name);
        end
    end
end
