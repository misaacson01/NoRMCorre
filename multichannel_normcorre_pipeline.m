
%set primary image filename (used to create a template during motion correction)
% imageName = 'C:\Users\misaa\Desktop\21-11-23 viral 3color satsuma and chameleon\21-11-23 chameleon 200depth tseries006.tif';
imageName = 'C:\Users\misaa\Desktop\21-11-23 viral 3color satsuma and chameleon\ds 21-11-23 satsuma 200depth002.tif';

%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [1 3]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels

%perform multichannel normcorre motion correction and create a new template
templateName = run_multichannel_normcorre('imageName',imageName,...
    'channelOptions',channel_options);

%%
%set additional filenames to register to the new template
imageNames{1} = 'C:\Users\misaa\Desktop\21-11-23 viral 3color satsuma and chameleon\ds 21-11-23 chameleon 200depth001.tif';

%change settings if desired
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels

%perform multichannel normcorre motion correction on each file
for i = 1:length(imageNames)
    run_multichannel_normcorre('imageName',imageNames{i},...
        'channelOptions',channel_options,'templateName',templateName);
end
