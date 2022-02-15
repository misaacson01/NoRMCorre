function imageName = folder2tif(folder)
%FUNCTION imageName = folder2tif(folder)
%
%combines a folder of .tif files into a single tif file

folder = 'C:\Users\misaa\Desktop\22-01-14 4color amazing data\Visual Cortex\dirtuning';
assert(isfolder(folder),'input must be a valid folder path')

files = dir(fullfile(folder,'*.tif'));
numFiles = length(files);

%generate a name for the full image file
[basePath, folderName] = fileparts(folder);
imageName = fullfile(basePath,[folderName '.tif']);

%if a file of that name already exist, append a timestamp
if exist(imageName,'file')==2
    imageName = fullfile(basePath,[folderName datestr(now,'yymmddHHMMSS') '.tif']);
end

opts_tiff.append = true;
opts_tiff.big = true;
opts_tiff.message = false;
        
prevstr = [];
for i = 1:numFiles
    images = read_file(fullfile(folder,files(i).name));
    saveastiff(images,imageName,opts_tiff);       
    str=[num2str(i), ' out of ', num2str(numFiles), ' files combined...'];
    refreshdisp(str, prevstr, i);
    prevstr=str; 
end
fprintf('done.\n')