function [status, templateName, tforms, corrcoefs] = run_multichannel_imregister(varargin)
% FUNCTION [status, templateName, tforms, corrcoefs] = run_multichannel_imregister(varargin)
%
% Processes multiple, multi-channel tif stacks through the normcorre motion
% correction scripts
%
% Optional function Inputs ('name':value pairs)
% 'imageName': full filename of .tif image series to be motion stabilized
% 'channel_options': 
%   channel_options.nch: number of channels in the source file
%   channel_options.chsh: channels to use for calculating shifts (e.g. [2 3 4])
%   channel_options.chrg: channels to register (e.g. [2 3 4])
%   channel_options.pr: projection type to use across channels ('max' or 'mean')
%   channel_options.save: whether to save the registered images, or just get the tforms (true or false)
%
% Output:
% status: text description of results (e.g. "success" or "bad template")
% templateName: directory/filename of the template generated from motion correction
% tforms: transformation matrices used to stabilize the frames
% corrcoefs: 
%
% example:
% \run_multichannel_imregister('sourceName','C:\image.tif','saveDir','C:\stabilized\','channel',2)

%get inputs
if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'imageName'
                imageName = varargin{i+1};
            case 'channelOptions'
                channel_options = varargin{i+1};
            case 'templateName'
                error('templateName parameter not yet supported')
            otherwise
                error('Input parameter name not recognized')
        end
    end
end

%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [1 2 3 4]; %channels to use for calculating shifts (e.g. [2 3 4])
channel_options.chrg = [1 2 3 4]; %channels to register (e.g. [2 3 4])
channel_options.chsv = [1 2 3 4]; %channels to save (e.g. [1 2 3 4])
channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')
channel_options.save = true; %whether to save the registered images, or just get the tforms (true or false)

%set input defaults
if ~exist('imageName','var')
    [imageFilename, imagePathname, ~] = uigetfile('*.tif','Select the .tif file to be motion stabilized');
    imageName = fullfile(imagePathname,imageFilename);
end
if ~exist('channel_options','var')
    channel_options.nch = 4; %number of channels in the source file
    channel_options.chsh = [1 2 3 4]; %channels to use for registering shifts
    channel_options.chrg = [1 2 3 4]; %channels to register
    channel_options.pr = 'max'; %projection type to use across channels
    channel_options.save = true; %whether to save the registered images or not (false just returns the last transform)
end
if ~isfield('channel_options','save')
    channel_options.save = true;
end

%check validity of inputs
[imagePathname,imageFilename,ext] = fileparts(imageName);
assert(~isempty(imagePathname),'imageName must include the full file path')
assert(strcmp(ext,'.tif'),'Source must be a .tif file')
assert(exist(imageName,'file')==2,['Cannot find file "' imageName '"']);

%additional settings
num_template_frames = 200;
opts_tiff.append = false;
opts_tiff.big = true;
opts_tiff.message = false;
[optimizer, metric] = imregconfig('multimodal');
gauss_filt = 0.5;

prevstr = '';
warning('off','images:regmex:registrationOutBoundsTermination')
warning('off','images:regmex:registrationFailedException');
warning('off','imageio:tiffmexutils:libtiffWarning');

%read file and get image size
data = read_file(imageName);
[h,w,num_total_frames] = size(data);
nch = channel_options.nch;
datatype = class(data);
num_frames = num_total_frames/nch;

data = reshape(data(:,:,1:num_frames*nch),[h w nch num_frames]);

%create projection of data
switch channel_options.pr
    case 'max'
        data_proj = max(data(:,:,channel_options.chsh,:),[],3);
    case 'mean'
        data_proj = mean(data(:,:,channel_options.chsh,:),[],3);
    otherwise
        error('projection type not supported')
end
data_proj = permute(data_proj,[1 2 4 3]);

%get mean of 1st N frames for initial template
if num_frames<200
    num_template_frames = num_frames;
end
template_proj = data_proj(:,:,1:num_template_frames);
template = mean(template_proj,3,'omitnan');

%register 1st N frames to create improve the initial template
for f = 1:num_template_frames
    str=['Registering first ' num2str(num_template_frames) ' frames to create an initial template... (' num2str(f) ')'];
    refreshdisp(str, prevstr, f);
    prevstr=str; 
    imr_tform = imregtform(template_proj(:,:,f),template,'affine',optimizer,metric);
    template_proj(:,:,f) = imwarp(template_proj(:,:,f),imr_tform,'OutputView',imref2d(size(template)));
end
template = mean(template_proj,3,'omitnan');
str=['Registering first ' num2str(num_template_frames) ' frames to create an initial template... done.\n'];
refreshdisp(str, prevstr, f); 
reg_tforms = nan(3,3,num_frames);


%% calculate shifts between all frames and the template
for f = 1:num_frames
    str=['Calculating ' num2str(num_frames) ' frame shifts... (' num2str(f) ')'];
    refreshdisp(str, prevstr, f);
    prevstr=str; 
    imr_tform = imregtform(imgaussfilt(data_proj(:,:,f),gauss_filt),template,'affine',optimizer,metric);
    reg_tforms(:,:,f) = imr_tform.T;
end
str = ['Calculating ' num2str(num_frames) ' frame shifts... done.\n'];
refreshdisp(str, prevstr, f);


%% optimize shifts
fprintf('Optimizing frame shifts... '); 
corrcoefs = nan(2,num_frames);
best_tmp_inds = nan(1,num_frames);
tforms = reg_tforms;
for f = 1:num_frames
    %get nearby registered shifts
    tmp_tforms = nan(3,3,7);
    tmp_cc = -1*ones(1,7);
    for r = -2:2
        frame_ind = f+r;
        tmp_ind = r+3;
        if frame_ind>=1 && frame_ind<=num_frames
            tmp_tforms(:,:,tmp_ind) = reg_tforms(:,:,frame_ind);
            imr_tform.T = reg_tforms(:,:,frame_ind); %warp based on tform
            tmp = imwarp(data_proj(:,:,f),imr_tform,'OutputView',imref2d(size(template)));
            tmp_cc(tmp_ind) = corr2(template,imgaussfilt(double(tmp),gauss_filt));
        end
    end
    
    %get average of best before/after shift
    tmp_ind = 6; 
    [~,best_before_tmp_ind] = max(tmp_cc(1:2));
    [~,best_after_tmp_ind] = max(tmp_cc(4:5));
    best_after_tmp_ind = best_after_tmp_ind + 3;
    tmp_tforms(:,:,tmp_ind) = mean(tmp_tforms(:,:,[best_before_tmp_ind best_after_tmp_ind]),3);
    if all(all(~isnan(tmp_tforms(:,:,tmp_ind))))
        imr_tform.T = tmp_tforms(:,:,tmp_ind); %warp based on tform
        tmp = imwarp(data_proj(:,:,f),imr_tform,'OutputView',imref2d(size(template)));
        tmp_cc(tmp_ind) = corr2(template,imgaussfilt(double(tmp),gauss_filt));
    end
    
    %no shift
    tmp_ind = 7;
    tmp_tforms(:,:,tmp_ind) = [1 0 0; 0 1 0; 0 0 1];
    tmp_cc(tmp_ind) = corr2(template,imgaussfilt(double(data_proj(:,:,f)),gauss_filt));
    corrcoefs(1,f) = tmp_cc(tmp_ind);
    
    %pick the best shift
    [~,best_tmp_ind] = max(tmp_cc);
    corrcoefs(2,f) = tmp_cc(best_tmp_ind);
    tforms(:,:,f) = tmp_tforms(:,:,best_tmp_ind);
    best_tmp_inds(f) = best_tmp_ind;
end
fprintf('done.\n');   


%% perform shifts
fprintf('Performing frame shifts... '); 
for f = 1:num_frames
    imr_tform.T = tforms(:,:,f);
    for c = 1:nch
        data(:,:,c,f) = imwarp(data(:,:,c,f),imr_tform,'OutputView',imref2d(size(template)));
    end
end
fprintf('done.\n');    


%% save registered images
if channel_options.save
    fprintf('Saving files... '); 
    savepath = fullfile(imagePathname,['IMR ' imageFilename]);
    if ~(exist(savepath,'dir')==7)
        mkdir(savepath);
    end
    %save template
    for c = channel_options.chsv
        saveastiff(cast(mean(data(:,:,c,:),4),datatype),fullfile(savepath,['AVG_CH' num2str(c) '_' imageFilename ext]),opts_tiff);
        saveastiff(cast(permute(data(:,:,c,:),[1 2 4 3]),datatype),fullfile(savepath,['CH' num2str(c) '_' imageFilename ext]),opts_tiff);
%         saveastiff(cast(max(data(:,:,c,:),[],4),datatype),fullfile(savepath,['MAX_CH' num2str(c) '_' imageFilename ext]),opts_tiff);
    end
    saveastiff(cast(template,datatype),fullfile(savepath,['TEMPLATE_' imageFilename ext]),opts_tiff);
    fprintf('done.\n');    
end

save(fullfile(savepath,['CorrCoef_before_after_' imageFilename '.mat']),'corrcoefs')
fprintf('Image registration complete.\n\n'); 
