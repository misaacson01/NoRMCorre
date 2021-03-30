clear all
%recommended: start by running each section one at a time

%% set file path of tif stack(s) or folder(s) of sequential tif images/stacks
tifNames{1} = 'C:\Users\Matthew\Documents\Schaffer-Nishimura Lab\Visual Stimulation\Data\21-03-26 100th viral gcamp spontaneous and directiontuning\directiontuning 180depth\'; 
tifNames{2} = 'C:\Users\Matthew\Documents\Schaffer-Nishimura Lab\Visual Stimulation\Data\21-03-26 100th viral gcamp spontaneous and directiontuning\21-03-26 100th viral gcamp spontaneous 180depth133.tif'; 
channels_to_save = 3; %set image channel (for multi-channel image series)


%% combine multiple tif files and folders into a single image file (if applicable)
[fullTifName, numTifFrames, fps] = combine_tifs(tifNames, channels_to_save);    


%% run normecorre function to save a motion-stabilized file
[savedir, ~, ~] = fileparts(fullTifName);
stabilizedFullTifName = run_normcorre('imageName',fullTifName,'saveDir',savedir);


%% add ROIs manually and extract fluorescence traces for each
draw_ROIs(stabilizedFullTifName, fps, numTifFrames); %run ROI-drawing GUI
