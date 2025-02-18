

%full filename of the multichannel .tif file to motion-correct
imageName = 'D:\SN Lab\Spinal Cord\data\051424_550615_4xzm_Lside_run4_spotA_00001.tif';

%settings for current image
channel_options.nch = 3; %number of channels in the source file
channel_options.chsh = [2]; %channels to use for calculating shifts (e.g. [2 3 4])
channel_options.chrg = [2]; %channels to register (e.g. [2 3 4])
channel_options.chsv = [2]; %channels to save (e.g. [1 2 3 4])
channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')

run_multichannel_imregister('imageName',imageName,'channelOptions',channel_options);