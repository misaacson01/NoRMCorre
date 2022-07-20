

%full filename of the multichannel .tif file to motion-correct
imageName = 'C:\Users\misaa\Desktop\07_15_22_Kailyn\071522_479576_256_spotE_airpuff2_014.tif';

%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for calculating shifts (e.g. [2 3 4])
channel_options.chrg = [2 3 4]; %channels to register (e.g. [2 3 4])
channel_options.chsv = [1 2 3 4]; %channels to save (e.g. [1 2 3 4])
channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')
channel_options.save = true; %whether to save the registered images, or just get the tforms (true or false)

run_multichannel_imregister('imageName',imageName,'channelOptions',channel_options);