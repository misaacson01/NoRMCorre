function savetiffframe(data,savepath,tagstruct,options)


%% Initialize options
assert(~isempty(data),'Image data is missing');
if nargin < 4 % Use default options
    options.append = false;
    options.overwrite = false;
end
if ~isfield(options, 'append')
    options.append = false; 
end
if ~isfield(options, 'overwrite')
    options.overwrite = false; 
end
if isfield(options, 'savechannels') && isfield(tagstruct, 'ImageDescription')
    tmp = tagstruct.ImageDescription;
    for i = 1:4
        ind = strfind(tmp,['savingChannel' num2str(i)]);
        if ~isempty(ind)
            tmp(ind+15) = num2str(options.savechannels(i));
        end
    end
    tagstruct.ImageDescription = tmp;
end


%% Overwrite check
if exist(savepath, 'file') && ~options.append && ~options.overwrite
    error('Save path/filename already exists, but neither append nor overwrite option was selected')
end


%% Save path configuration
[pathstr, fname, fext] = fileparts(savepath);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
end

%% Write image data to a file
file_opening_error_count = 0;
while ~exist('tfile', 'var')
    try
        if ~options.append % Make a new file
            tfile = Tiff(savepath, 'w');
        else
            if ~exist(savepath, 'file') % Make a new file
                tfile = Tiff(savepath, 'w');
            else % Append to an existing file
                tfile = Tiff(savepath, 'r+');
                while ~tfile.lastDirectory(); % Append a new image to the last directory of an exiting file
                    tfile.nextDirectory();
                end
                tfile.writeDirectory();
            end
        end
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                error('Could not open tiff file')
            end
        end
    end
end

tfile.setTag(tagstruct);
tfile.write(data);
tfile.close();


end