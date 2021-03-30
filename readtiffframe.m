function [data, tagstruct, channels] = readtiffframe(filename, frameNum)
% FUNCTION [data, tagstruct, channels] = readtiffframe(path, frameNum)
%
% reads a single frame of a Tiff file including all metadata commonly used
% by scanimage
%
% INPUTS
%filename: filename of .tif file
%frameNum: frame number of tif file to read (for multi-image Tiff files)
%
%OUTPUTS
%data: image data from selected frame of Tiff file
%tagstruct: metadata associated with the selected frame
%channels: channels [1-4] that are being saved in the entire Tiff files


%% Read specified frame of tiff file
%turn unnecessary Tiff library warnings off
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
warning('off','MATLAB:imageio:tiffmexutils:libtiffWarning')
t = Tiff(filename,'r');
setDirectory(t,frameNum);
data = imread(filename,frameNum);


%% Get all image tags and save in tagstruct (i.e. image metadata)
tagstruct.SubFileType = getTag(t,'SubFileType');
tagstruct.Photometric = getTag(t,'Photometric');
tagstruct.ImageLength = getTag(t,'ImageLength');
tagstruct.ImageWidth = getTag(t,'ImageWidth');
tagstruct.RowsPerStrip = getTag(t,'RowsPerStrip');
tagstruct.BitsPerSample = getTag(t,'BitsPerSample');
tagstruct.Compression = getTag(t,'Compression');
tagstruct.SampleFormat = getTag(t,'SampleFormat');
tagstruct.SamplesPerPixel = getTag(t,'SamplesPerPixel');
tagstruct.PlanarConfiguration = getTag(t,'PlanarConfiguration');
tagstruct.Orientation = getTag(t,'Orientation');

%try to get ImageDescription tag (contains most of scanimage metadata)
channels = nan([1,4]);
try
    tagstruct.ImageDescription = getTag(t,'ImageDescription');
    tmp = tagstruct.ImageDescription;
    for i = 1:4
        ind = strfind(tmp,['savingChannel' num2str(i)]);
        if ~isempty(ind)
            channels(i) = str2double(tmp(ind+15));
        else
            error(['Cannot find savingChannel' num2str(i) ' in ImageDescription'])
        end
    end
catch
    error('Cannot retrieve ImageDescription from tiff header')
end

%% close Tiff file
t.close();
%turn warnings back on
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')
warning('on','MATLAB:imageio:tiffmexutils:libtiffWarning')

end