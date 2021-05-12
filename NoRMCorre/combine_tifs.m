function [combined_filename, num_frames, fps] = combine_tifs(tifNames, channels_to_save)
% FUNCTION combined_filename = combine_tifs(tifNames, channels_to_save)
%
% combines all .tif files (and folders of files) into a single .tif image 
% stack
% 
%INPUTS
%tifNames: cell array of .tif filenames and/or directories containing .tif 
%          files to be combined
%channels_to_save: (optional) image channels to save in combined file
%                  (if left blank, all channels will be saved)
%
%OUTPUTS:
%combined_filename: full filename of combined .tif file
%
%example:
%combined_filename = combine_tifs('C:/images/',3); %save channel 3 only of 
%                    all tif files in images folder
                    
                    
%% check that inputs make sense (and get save directory for combined file)
%check that all files/folders exist
numTifNames = length(tifNames);
for n = 1:numTifNames
    assert(exist(tifNames{n},'dir')==7||exist(tifNames{n},'file')==2,['Cannot find "' tifNames{n} '"']);
end

%set a filename for the combined file
[combined_dir, ~, ~] = fileparts(tifNames{1});
combined_filename = fullfile(combined_dir,'combined_images.tif'); %name for combined files
num_frames = nan(numTifNames,1);
fps = nan(numTifNames,1);
assert(~(exist(combined_filename,'file')==2),'combined filename already exists');

%check/set optional 2nd input
if nargin<2
    channels_to_save = [1 2 3 4];
end
assert(channels_to_save>0 & channels_to_save<5 & all(mod(channels_to_save,1)==0),...
    'channels_to_save must be integers from 1-4')
channels_to_save = sort(channels_to_save); %make sure it's in ascending order
options.savechannels = [0 0 0 0]; %convert save channels to binary array
options.savechannels(channels_to_save) = 1;
if length(channels_to_save)==1
    savestr = ['channel ' num2str(channels_to_save)];
else
    savestr = ['channels [' num2str(channels_to_save) ']' ];
end
options.append = false; 


%% loop for all filenames/folder to combine all .tif files
for n = 1:numTifNames
    if exist(tifNames{n},'dir')==7 % for .tif file folders
        %find all tif files and sort by date created
        folderTifNames = dir([tifNames{n} '\*.tif']);
        [~,idx] = sort([folderTifNames.datenum]);
        folderTifNames = folderTifNames(idx);
        num_tifs = length(folderTifNames);
        if num_tifs>size(num_frames,2)
            num_frames(:,end+1:num_tifs) = nan;
        end
        
        %load every image frame and save into a new, single tiff stack
        fprintf('%s',['Folder: '  tifNames{n}]);
        fprintf('\n');
        fprintf(['Combining ' savestr ' from ' num2str(num_tifs) ' tif files... (1)']);
        backspaces = 2;
        for t = 1:num_tifs
            fprintf(1,[repmat('\b',[1 backspaces]) num2str(t) ')']);
            backspaces = 1 + numel(num2str(t));

            curName = folderTifNames(t).name;
            curFullName = fullfile(folderTifNames(t).folder,curName);
            if t==1
                curTiff = Tiff(curFullName,'r');
                fps(n) = read_tiffstate(curTiff,'frameRate');
            end
            
            info = imfinfo(curFullName);
            num_frames(n,t) = length(info);
            [~, ~, channels_bin] = readtiffframe(curFullName, 1);
            channels = find(channels_bin==1);
            for f = 1:num_frames(n,t)
                cur_channel = channels(mod(f-1,length(channels))+1);
                if any(channels_to_save==cur_channel)
                    [data, tagstruct, ~] = readtiffframe(curFullName, f);
                    savetiffframe(data,combined_filename,tagstruct,options);
                    if options.append == false
                        options.append = true; %after first saved frame, append the rest
                    end
                end
            end
            num_frames(n,t) = length(channels_to_save)*num_frames(n,t)/length(channels);
        end
        fprintf(1,[repmat('\b',[1 1+backspaces]) 'done.\n']);
        
    else %for individual .tif files

        %print names of files to be combined
        fprintf('%s',['File: '  tifNames{n}]);
        fprintf('\n');
        fprintf(['Combining ' savestr '... ']);
        
        curFullName = tifNames{n};
        curTiff = Tiff(curFullName,'r');
        fps(n) = read_tiffstate(curTiff,'frameRate');
            
        info = imfinfo(curFullName);
        num_frames(n,1) = length(info);
        [~, ~, channels_bin] = readtiffframe(curFullName, 1);
        channels = find(channels_bin==1);
        for f = 1:num_frames(n,1)
            cur_channel = channels(mod(f-1,length(channels))+1);
            if any(channels_to_save==cur_channel)
                [data, tagstruct, ~] = readtiffframe(curFullName, f);
                savetiffframe(data,combined_filename,tagstruct,options);
                if options.append == false
                    options.append = true; %after first saved frame, append the rest
                end
            end
        end
        num_frames(n,1) = length(channels_to_save)*num_frames(n,1)/length(channels);
        fprintf('done.\n');
    end
    
end %end tifName loop
   
end %end function