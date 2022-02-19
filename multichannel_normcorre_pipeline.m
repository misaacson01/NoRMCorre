%% register a multichannel image series, creating a new template
%set primary image filename (used to create a template during motion correction)
% imageName = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\920_00001x.tif';
imageName = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\1030_00001.tif';


%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels

%perform multichannel normcorre motion correction and create a new template
[status, templateName, ~] = run_multichannel_normcorre('imageName',imageName,...
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
    [status, ~, ~] = run_multichannel_normcorre('imageName',imageNames{i},...
        'channelOptions',channel_options,'templateName',templateName);
end


%% Take 2 individually-registered datasets and register them together using their templates
imageFolder1 = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\MC 920_00001';
template1 = 'TEMPLATE_920_00001';
imageFolder2 = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\MC dirtuning';
template2 = 'TEMPLATE_dirtuning';

%register template2 to template1 with imregtform (better for big motion)
template1Name = fullfile(imageFolder1,[template1 '.tif']);
template2Name = fullfile(imageFolder2,[template2 '.tif']);
[optimizer, metric] = imregconfig('multimodal');
fixed = imread(template1Name);
moving = imread(template2Name);
imr_tform = imregtform(moving,fixed,'affine',optimizer,metric);
movingRegistered = imwarp(moving,imr_tform,'OutputView',imref2d(size(fixed)));

%get correlation between registered template2 and template1
kept_mask = imwarp(ones(size(fixed)),imr_tform,'OutputView',imref2d(size(fixed)));
fixed_masked = fixed; 
fixed_masked(~kept_mask) = 0;
R = corr2(fixed,moving);
R2 = corr2(fixed_masked,movingRegistered);

%save registered template2
template2Name_reg = fullfile(imageFolder1,[template2 '_reg2_' template1 '.tif']);
opts_tiff.append = true;
opts_tiff.big = true;
opts_tiff.message = false;
saveastiff(movingRegistered,template2Name_reg,opts_tiff);  

%register template2 a 2nd time using normcorre (better for small, in-frame motion)
channel_options.nch = 1;
channel_options.chsh = 1;
channel_options.pr = 'max';
channel_options.save = false;
[status, ~, nrm_tform] = run_multichannel_normcorre('imageName',template2Name_reg,...
        'channelOptions',channel_options,'templateName',template1Name);
movingRegistered2 = imwarp(movingRegistered,nrm_tform,'cubic','FillValues',0);
kept_mask = imwarp(kept_mask,nrm_tform,'cubic','FillValues',0);
fixed_masked2 = fixed; 
fixed_masked2(~kept_mask) = 0;
R3 = corr2(fixed_masked,movingRegistered2);
figure()
subplot(2,2,1)
image(fixed_masked)
subplot(2,2,2)
image(fixed_masked2)
subplot(2,2,3)
image(movingRegistered)
subplot(2,2,4)
image(movingRegistered2)

%as user to continue if template2 could be well aligned to template1
if R2>R3
    answer = input(['Correlation improved from R=' num2str(R) ' to R=' num2str(R2) ' using imregtform only, register the entire folder this way? (y/n): '],'s');
else
    answer = input(['Correlation improved from R=' num2str(R) ' to R=' num2str(R3) ' using both imregtform and normcorre, register the entire folder this way? (y/n): '],'s');
end
if strncmpi(answer,'y',1)
    %warp every image file in folder2
    files2 = dir(fullfile(imageFolder2,'*.tif'));
    numFiles = length(files2);
    prevstr = [];
    saveFolder = fullfile(imageFolder2,['registered to ' template1]);
    if isfolder(saveFolder)
        error(['"' saveFolder '" folder already exists'])
    else
        mkdir(saveFolder);
    end
    for f = 1:numFiles
        images = read_file(fullfile(imageFolder2,files2(f).name));
        for i = 1:size(images,3)
            images(:,:,i) = imwarp(images(:,:,i),imr_tform,'OutputView',imref2d(size(fixed)));
            if R3>R2
                images(:,:,i) = imwarp(images(:,:,i),nrm_tform,'cubic','FillValues',0);
            end
        end
        saveastiff(images,fullfile(saveFolder,files2(f).name),opts_tiff);       
        str=[num2str(f), ' out of ', num2str(numFiles), ' files warped to template1...'];
        refreshdisp(str, prevstr, f);
        prevstr=str; 
    end
end
fprintf('Done.\n')