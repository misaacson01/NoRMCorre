
clear all
%recommended: start by running each section one at a time


%% set file path of tif stack(s) or folder(s) of sequential tif images/stacks
dir = 'C:\Users\Matthew\Documents\Schaffer-Nishimura Lab\Visual Stimulation\Data\21-08-17 wt 3color mira and satsuma';
tifNames{1} = fullfile(dir,'21-08-17 wt 3color mira 910nm 80fr004.tif'); 
tifNames{2} = fullfile(dir,'21-08-17 wt 3color satsuma 1030nm 80fr005.tif'); 
channels_to_save = [1 3]; %set image channels to save
channel_to_register = 3; %choose best image channel for registrations
channels_to_apply = 1; %choose all non-registration channels to apply motion correction to


%% combine multiple tif files and folders into a single image file for each channel to save (if applicable)
[fullTifNames, numTifFrames, fps] = combine_tifs(tifNames, channels_to_save);    


%% run normecorre function to register one file
reg_idx = channels_to_save==channel_to_register;
[savedir, ~, ~] = fileparts(fullTifNames{1});
[stabilizedFullTifName, shifts, col_shift, options] = run_normcorre('imageName',fullTifNames{reg_idx},'saveDir',savedir,'numChannels',1);
fprintf('Finished registering.\n')


%% apply shifts from registered file to all others
fprintf('Applying shifts...\n')
stabilizedfullTifNames = cell(length([channel_to_register channels_to_apply]),1);
stabilizedfullTifNames{1} = stabilizedFullTifName;
num_to_apply = length(channels_to_apply);
for s = 1:num_to_apply
    apply_idx = channels_to_save==channels_to_apply(s);
    options.tiff_filename = fullTifNames{apply_idx};
    options.output_filename = fullfile(dir,['normcorre_combined_images_ch' num2str(channels_to_apply(s)) '.tif']);
    M_final = apply_shifts(fullTifNames{apply_idx},shifts,options,0,0,0,col_shift);
    stabilizedfullTifNames{1+s} = M_final;
end
fprintf('Finished applying shifts.\n')


%% add ROIs manually and extract fluorescence traces for each
draw_ROIs(stabilizedFullTifName, fps, numTifFrames); %run ROI-drawing GUI

