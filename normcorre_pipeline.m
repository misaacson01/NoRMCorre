clear all
%recommended: start by running each section one at a time

%% set file path of tif stack(s) or folder(s) of sequential tif images/stacks
% tifNames{1} = 'C:\Users\Matthew\Documents\Schaffer-Nishimura Lab\Visual Stimulation\Data\21-04-29 viral spatial\spatialtuning\'; 
tifNames{1} = 'D:\SN Lab\Gut\Imaging\C2-Zeiss 20X zoom1-2 256X256_00001.tif';
channels_to_save = 2; %3;%set image channel (for multi-channel image series)


%% combine multiple tif files and folders into a single image file (if applicable)
[fullTifName, numTifFrames, fps] = combine_tifs(tifNames, channels_to_save);    


%% run normecorre function to save a motion-stabilized file
fullTifName = 'D:\SN Lab\Gut\Imaging\dataset3\C2-20X zoom 1 cap treat video3_00001.tif';
[savedir, ~, ~] = fileparts(fullTifName);
stabilizedFullTifName = run_normcorre('imageName',fullTifName,'saveDir',savedir);


%% add ROIs manually and extract fluorescence traces for each
draw_ROIs(stabilizedFullTifName, fps, numTifFrames); %run ROI-drawing GUI
