
file_path = "/Users/ryan.gorzek/Library/CloudStorage/GoogleDrive-ryan.gorzek@gmail.com/Shared drives/Astrocytes - Protocadherins/";
addpath(genpath(fullfile(file_path)));

subj_ids = ["john01", "john02", "john03"];
plane_ids = ["000", "001", "002"];
eye_ids = ["000", "001"]; % contra, ipsi
for subj = subj_ids
    for plane = plane_ids
        exp_id = char(strcat(subj, "_", plane));
        fprintf("Processing %s...", exp_id);
        exp_path = char(strcat(file_path, filesep, subj, filesep, exp_id, filesep));
        plane_path = strcat(exp_path, "suite2p", filesep, "plane0", filesep);
        % Generate .align, .segment, and .signals1 files. 
        sbxsuite2sbx(strcat(plane_path, "Fall"), strcat(exp_path, exp_id));
        % Generate .signals file.
        sbxf2spks(strcat(exp_path, exp_id));
        % Split .signals file for the two eyes.
        sbxsplitsuite(strcat(exp_path, exp_id));
        % Calculate tuning results and generate .orisf file for each eye.
        for eyeball = eye_ids
            eye_id = char(strcat(exp_path, exp_id, "_", eyeball, "_suite"));
            sbxorisf_new(eye_id);
        end
    end
end
