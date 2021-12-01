function templateName = run_multichannel_normcorre(varargin)
% FUNCTION saveName = run_multichannel_normcorre(varargin)
%
% Processes multiple, multi-channel tif stacks through the normcorre motion
% correction scripts
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
            case 'channelOptions'
                channel_options = varargin{i+1};
            case 'templateName'
                templateName = varargin{i+1};
            otherwise
                error('Input parameter name not recognized')
        end
    end
end
if ~exist('imageName','var')
    [imageFilename, imagePathname, ~] = uigetfile('*.tif','Select the .tif file to be motion stabilized');
    imageName = fullfile(imagePathname,imageFilename);
end
if ~exist('channel_options','var')
    channel_options.nch = 4; %number of channels in the source file
    channel_options.chsh = [1 2 3 4]; %channels to use for registering shifts
    channel_options.pr = 'max'; %projection type to use across channels
end
[imagePathname,imageFilename,ext] = fileparts(imageName);
assert(~isempty(imagePathname),'imageName must include the full file path')
assert(strcmp(ext,'.tif'),'Source must be a .tif file')
assert(exist(imageName,'file')==2,['Cannot find file "' imageName '"']);

info = imfinfo(imageName);
image_size = [info(1).Height, info(1).Height];

%apply motion correction using normCorre algorithm
normcorre_options = NoRMCorreSetParms('d1',image_size(1),'d2',image_size(2),'grid_size',options.grid_size,'init_batch',options.init_batch,...
                'overlap_pre',options.overlap_pre,'mot_uf',options.mot_uf,'bin_width',options.bin_width,'max_shift',options.max_shift,...
                'max_dev',options.max_dev,'us_fac',options.us_fac,'output_type','tif','tiff_filename',imageName);
if exist('templateName','var') && ischar(templateName) && length(templateName)>4
    template = imread(templateName);
    normcorre_options.upd_template = false; %don't change the template
    [~,~,~,~,~] = multichannel_normcorre(imageName,normcorre_options,channel_options,template);
    templateName = '';
else
    normcorre_options.upd_template = true; %create a template
    [~,~,template,~,~] = multichannel_normcorre(imageName,normcorre_options,channel_options);
    %save template
    templateName = fullfile(imagePathname,['MC ' imageFilename],['TEMPLATE_' imageFilename ext]);
    imwrite(uint16(template),templateName);
end

end