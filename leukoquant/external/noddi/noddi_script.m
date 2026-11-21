function matlab_script(dwi_path, brain_mask_path, bvecs_path, bvals_path, output_path, cores)    

    cd(fullfile(output_path));

    data_filename_hdr = strrep(dwi_path, '.nii.gz', '.hdr');
    data_filename_img = strrep(dwi_path, '.nii.gz', '.img');

    brain_mask_filename_img = strrep(brain_mask_path, '.nii.gz', '.img');
    brain_mask_filename_hdr = strrep(brain_mask_path, '.nii.gz', '.hdr');
        
    convertNiiToHdr(dwi_path);
    convertNiiToHdr(brain_mask_path);

    bvals_length_fix(bvals_path, bvecs_path);

    noddi_bvecs_path = fullfile(output_path, "noddi.bvecs");
    noddi_bvals_path = fullfile(output_path, "noddi.bvals");

    transpose_file(bvecs_path, noddi_bvecs_path);
    transpose_file(bvals_path, noddi_bvals_path);
    
    noddi_roi_path = fullfile(output_path, 'NODDI_roi.mat');
    CreateROI(data_filename_hdr, brain_mask_filename_hdr, noddi_roi_path);

    protocol = FSL2Protocol(noddi_bvals_path, noddi_bvecs_path);
    noddi = MakeModel('WatsonSHStickTortIsoV_B0');
    
    if nargin < 6
        cores = 1;
    end
    
    noddi_roi_filename = 'NODDI_roi.mat';
    fitted_params_filename = 'FittedParams.mat';

    %batch_fitting_single(noddi_roi_filename, protocol, noddi, 'FittedParams.mat');    
    batch_fitting(noddi_roi_filename, protocol, noddi, fitted_params_filename, cores);
    
    output_prefix = '';
    SaveParamsAsNIfTI(fitted_params_filename, noddi_roi_filename, brain_mask_filename_hdr, output_prefix);
    
    delete_file(dwi_path);
    delete_file(data_filename_img);
    delete_file(data_filename_hdr);
    delete_file(brain_mask_filename);
    delete_file(brain_mask_filename_hdr);
    delete_file(brain_mask_filename_img);
end

function convertNiiToHdr(niiFile)
    % Get the data and header information
    data = niftiread(niiFile);
    header = niftiinfo(niiFile);

    % Remove the file extension from the input filename
    [~, basename, ~] = fileparts(niiFile);

    % Save the header information to the Analyze header file
    niftiwrite(data,basename,header, 'Combined', false);
    
    delete_file(niiFile);
    
    disp(['Conversion completed for ' niiFile]);
end

% Create a function that counts how many lines are in bvecs and bvals and repeats the contents of bvals to match the number of lines in bvecs
function bvals_length_fix(bvals, bvecs)
    % Read the bvals and bvecs files
    bvals = load(bvals);
    bvecs = load(bvecs);
    
    % Get the number of lines in bvecs
    bvecs_length = size(bvecs, 1);
    
    % Get the number of lines in bvals
    bvals_length = size(bvals, 1);
    
    % If the number of lines in bvals is less than the number of lines in bvecs, repeat the contents of bvals to match the number of lines in bvecs
    if bvals_length < bvecs_length
        bvals = repmat(bvals, bvecs_length/bvals_length, 1);
    end
    
    % Save the fixed bvals to the same file
    save(bvals, 'bvals', '-ascii');
end

function transpose_file(file_temp, new_filename)
    data_temp = transpose(load(file_temp));
    save(new_filename, 'data_temp', "-ascii");
end
    
function delete_file(filename)
    delete(filename);
    disp(['Deleting file ' filename]);  
end
    

