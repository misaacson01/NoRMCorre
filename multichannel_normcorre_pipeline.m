%% register a multichannel image series, creating a new template
%set primary image filename (used to create a template during motion correction)
% imageName = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\920_00001x.tif';
imageName = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\1030_00001.tif';


%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels

%perform multichannel normcorre motion correction and create a new template
[status, templateName] = run_multichannel_normcorre('imageName',imageName,...
    'channelOptions',channel_options);


%% register additional multichannel image series to the generated template
%set additional filenames to register to the new template
clear imageNames
% imageNames{1} = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\1030_00001.tif';
imageNames{1} = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\dirtuning.tif';
templateName = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\MC 920_00001\TEMPLATE_920_00001.tif';

%change settings if desired
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels

%perform multichannel normcorre motion correction on each file
for i = 1:length(imageNames)
    [status, ~] = run_multichannel_normcorre('imageName',imageNames{i},...
        'channelOptions',channel_options,'templateName',templateName);
end


%% Take 2 individually-registered datasets and register them together using their templates
imageFolder1 = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\MC 920_00001';
template1 = 'TEMPLATE_920_00001';
imageFolder2 = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\MC dirtuning';
template2 = 'TEMPLATE_dirtuning';

[optimizer, metric] = imregconfig('multimodal');
fixed = imread(fullfile(imageFolder1,[template1 '.tif']));
moving = imread(fullfile(imageFolder2,[template2 '.tif']));
tform = imregtform(moving,fixed,'affine',optimizer,metric);
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
kept_mask = imwarp(ones(size(fixed)),tform,'OutputView',imref2d(size(fixed)));
fixed_masked = fixed; fixed_masked(~kept_mask) = 0;
R = corr2(fixed,moving);
R2 = corr2(fixed_masked,movingRegistered);
answer = input(['Correlation improved from R=' num2str(R) ' to R=' num2str(R2) ', register the entire folder? (y/n): '],'s');
if strncmpi(answer,'y',1)
    %warp every image file in folder2
    files2 = dir(fullfile(imageFolder2,'*.tif'));
    numFiles = length(files2);
    prevstr = [];
    opts_tiff.append = true;
    opts_tiff.big = true;
    opts_tiff.message = false;
    saveFolder = fullfile(imageFolder2,['registered to ' template1]);
    if isfolder(saveFolder)
        error(['"' saveFolder '" folder already exists'])
    else
        mkdir(saveFolder);
    end
    for f = 1:numFiles
        images = read_file(fullfile(imageFolder2,files2(f).name));
        for i = 1:size(images,3)
            images(:,:,i) = imwarp(images(:,:,i),tform,'OutputView',imref2d(size(fixed)));
        end
        saveastiff(images,fullfile(saveFolder,files2(f).name),opts_tiff);       
        str=[num2str(f), ' out of ', num2str(numFiles), ' files warped to template1...'];
        refreshdisp(str, prevstr, f);
        prevstr=str; 
    end
end
fprintf('Done.\n')