clear all
%recommended: start by running each section one at a time

%% set file path of tif stack or folder of sequential tif images/stacks
tifName = 'C:\image.tif'; %for example, a 4-channel tif stack
channels_to_save = 3; %set image channel (for multi-channel image series)


%% combine multiple tif files (if tifName is a folder rather than a single image file)
if isfolder(tifName) %combine tif images in folder for motion correction
    tifName = combine_tifs(tifName, channels_to_save);    
end


%% run normecorre function to save a motion-stabilized file
[savedir, ~, ~] = fileparts(tifName);
t = Tiff(tifName,'r');
fps = read_tiffstate(t,'frameRate');
stabilized_filename = run_normcorre('imageName',tifName,'saveDir',savedir);


%% add ROIs manually and extract fluorescence traces for each
draw_ROIs(stabilized_filename, fps); %run ROI-drawing GUI
