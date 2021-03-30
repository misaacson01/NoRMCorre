function saveName = run_normcorre(varargin)
% FUNCTION saveName = run_normcorre(varargin)
%
% Reads a .tif file and saves a motion-stabilized version
%
% Optional function Inputs ('name':value pairs)
% 'imageName': full filename of .tif image series to be motion stabilized
% 'saveDir': directory to save stabilized video
% 'channel': image channel of source video to stabilize (for mult-channel stacks)
%
% Output:
% saveName: directory/filename of stabilized video
%
% example:
% stabilized_filename = run_normcorre('sourceName','C:\image.tif','saveDir','C:\stabilized\','channel',2)

%options for normcorre algorithm
options.grid_size = [32,32];
options.init_batch = 200;
options.overlap_pre = 8;
options.mot_uf = 4;
options.bin_width = 200;
options.max_shift = 24;
options.max_dev = 8;
options.us_fac = 50;
%get inputs
if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'imageName'
                imageName = varargin{i+1};
            case 'saveDir'
                saveDir = varargin{i+1};
            case 'channel'
                channel = varargin{i+1};
            otherwise
                error('Input parameter name not recognized')
        end
    end
end
if ~exist('imageName','var')
    [imageFilename, imagePathname, ~] = uigetfile('*.tif','Select the .tif file to be motion stabilized');
    imageName = fullfile(imagePathname,imageFilename);
end
if ~exist('saveDir','var')
    saveDir = uigetdir('Choose folder to save stabilized .tif file');
end

%turn off unnecessary Tiff library warning
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning'); 
warning('off','MATLAB:imageio:tiffmexutils:libtiffWarning'); 

[~,imageFilename,ext] = fileparts(imageName);
assert(exist(imageName,'file')==2,['Cannot find file "' imageName '"']);
assert(exist(saveDir,'dir')==7,['Cannot find directory "' saveDir '"']);
assert(strcmp(ext,'.tif'),'Source must be a .tif file')

%get image info using Matlab Tiff library
t = Tiff(imageName,'r');
image_size = [getTag(t,'ImageLength') getTag(t,'ImageWidth')];
saved_channels = [read_tiffstate(t,'savingChannel1') read_tiffstate(t,'savingChannel2') ...
    read_tiffstate(t,'savingChannel3') read_tiffstate(t,'savingChannel4')];
num_channels = nansum(saved_channels);
if num_channels==0 %if savedchannels cannot be read, read all frames
    num_channels = 1; 
end
info = imfinfo(imageName);
num_images = length(info);
num_frames = num_images/num_channels;
assert(mod(num_frames,1)==0,'Incorrect number of channels and/or frames');

image_series = nan([image_size, num_frames]);
saveName = fullfile(saveDir,['normcorre_' imageFilename ext]);

%set channel (if multi-channel image stack) to stabilize
if ~exist('channel','var')
    if num_channels==1
        channel = 1;
    else
        channel = input(['Which channel (1-' num2str(num_channels) ') do you want to stabilize? ']);
    end
end
assert(any(channel==[1 2 3 4]),['Channel must be 1-' num2str(num_channels)]);

%load image data
chan_inds = channel:num_channels:num_images;
if num_channels==1
    fprintf(1,['Loading ' num2str(num_frames) ' images for motion correction. (1)']);
else
    fprintf(1,['Loading ' num2str(num_frames) ' images (from channel ' num2str(channel) ...
    ' of ' num2str(num_channels) '-channel stack) for motion correction. (1)']);
end
backspaces = 2;
for f = 1:num_frames
    if mod(f,10)==0
        fprintf(1,[repmat('\b',[1 backspaces]) num2str(f) ')']);
        backspaces = 1 + numel(num2str(f));
    end
    image_series(:,:,f) = imread(imageName,chan_inds(f));
end
fprintf([repmat('\b',[1 backspaces+1]) 'done.\n']);

%apply motion correction using normCorre algorithm
normcorre_options = NoRMCorreSetParms('d1',image_size(1),'d2',image_size(2),'grid_size',options.grid_size,'init_batch',options.init_batch,...
                'overlap_pre',options.overlap_pre,'mot_uf',options.mot_uf,'bin_width',options.bin_width,'max_shift',options.max_shift,...
                'max_dev',options.max_dev,'us_fac',options.us_fac,'output_type','tif','tiff_filename',saveName);
[~,~,~,~,~] = normcorre(image_series,normcorre_options);

%turn warning back on
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('on','MATLAB:imageio:tiffmexutils:libtiffWarning');

end