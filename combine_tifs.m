function combined_filename = combine_tifs(folder, channels_to_save)
% FUNCTION combined_filename = combine_tifs(folder, channels_to_save)
%
% combines all .tif files in a folder into a single .tif image stack.
% 
%INPUTS
%folder: directory containing .tif files to be combined
%channels_to_save: (optional) image channels to save in combined file
%(if left blank, all channels will be saved)
%
%OUTPUTS:
%combined_filename: full filename of combined .tif file
%
%example:
%combined_filename = combine_tifs('C:/images/',3); %save channel 3 only


%check that inputs make sense
assert(exist(folder,'dir')==7,['Cannot find folder "' folder ]);
if nargin<2
    channels_to_save = [1 2 3 4];
end
assert(channels_to_save>0 & channels_to_save<5 & all(mod(channels_to_save,1)==0),...
    'channels_to_save must be integers from 1-4')
channels_to_save = sort(channels_to_save); %make sure it's in ascending order

%set a filename for the combined file
combined_filename = fullfile(folder,'combined_images.tif'); %name for combined files

%find all tif files and sort by date created
tifNames = dir([folder '\*.tif']);
[~,idx] = sort([tifNames.datenum]);
tifNames = tifNames(idx);
num_tifs = length(tifNames);
num_frames = nan(1,num_tifs);

%check that files to be combined does not already include a combined file
for f = 1:num_tifs
    if strcmp(tifNames(f).name,'combined_images.tif')
        error('folder already includes a "combined_images.tif" file. Exiting script.')
    end
end

%print names of files to be combined
fprintf('Files to be combined:\n');
for f = 1:num_tifs
    fprintf([tifNames(f).name '\n']);
end
fprintf('\n');

%set options for saving combined file
options.savechannels = [0 0 0 0]; %convert save channels to binary array
options.savechannels(channels_to_save) = 1;
if length(channels_to_save==1)
    savestr = ['channel ' num2str(channels_to_save)];
else
    savestr = ['channels [' num2str(channels_to_save) ']' ];
end
options.append = false; 

%load every image frame and save into a new, single tiff stack
fprintf(['Combining ' savestr ' from ' num2str(num_tifs) ' tif files... (1)']);
backspaces = 2;
for t = 1:num_tifs
    fprintf(1,[repmat('\b',[1 backspaces]) num2str(t) ')']);
    backspaces = 1 + numel(num2str(t));

    curName = tifNames(t).name;
    curFullName = fullfile(folder,curName);
    info = imfinfo(curFullName);
    num_frames(t) = length(info);
    [~, ~, channels_bin] = readtiffframe(curFullName, 1);
    channels = find(channels_bin==1);
    for f = 1:num_frames(t)
        cur_channel = channels(mod(f-1,length(channels))+1);
        if any(channels_to_save==cur_channel)
            [data, tagstruct, ~] = readtiffframe(curFullName, f);
            savetiffframe(data,combined_filename,tagstruct,options);
            if options.append == false
                options.append = true; %after first saved frame, append the rest
            end
        end
    end
    num_frames(t) = length(channels_to_save)*num_frames(t)/length(channels);
end
fprintf(1,[repmat('\b',[1 1+backspaces]) 'done.\n']);
    
end