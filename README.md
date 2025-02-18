# Notes on this fork

This is a fork of the [NoRMCorre repo](https://github.com/flatironinstitute/NoRMCorre) developed by the Flatiron Institute, intended to more easily fit into S-N Lab image processing pipelines and to add some additional features to the motion correction tool. The [draw_ROIs repo](https://github.com/sn-lab/draw_ROIs) may be required to work with some of these scripts.

### Basic .tif file registration with default settings
To apply the basic NoRMCorre motion correction process with default settings to a scanimage .tif file, first identify whether the .tif file is single- or multi-channel. If it is single-channel, simply call the `run_normcorre` script from MATLAB with the full image filename and a save filename, like this:
`run_normcorre('imageName','c:\unstabilized.tif','saveDir','c:\stabilized.tif')`
If the .tif file is multichannel, include additional inputs to the function to specify the number of channels in the file and the channel you would like to correct, like this:
`run_normcorre('imageName','c:\unstabilized.tif','saveDir','c:\stabilized.tif','channel',2,'num_channels',3)`

### Registering multichannel .tif files
In addition to the basic NoRMCorre function which stabilizes a single imaging channel, this repository also includes a `run_multichannel_normcorre` function which can use one or multiple channels to stabilize all channels in a multi-channel .tif file. This is accomplished by taking a projection of all channels that will be used to stabilize the image, creating a stable template based on that projection, and registering the multi-channel projection of each frame of the .tif file to that template. After shifts have been identified which best stabilize the projection to the template, those shifts are then individually applied to one or a number of channels in the multi-channel .tif file, so that all are stabilized similarly. This process can take advantage of the fact that some channels tend to prodce higher quality images, whereas other channels may can be more sparse or noisy which makes them difficult to stabilize. By apply the shifts gained from registering high-quality channels and applying those shifts to all channels, noisy channels can be registered just as well as the best channels. 

To perform multi-channel image registration, check the examples located in the `multichannel_normcorre_pipeline.m` script.

### Manually selecting stable frames to improve the initial template

To perform motion correction of an image series, first you need a "template" to register each frame of the image series onto. A template could be generated simply by registering a subset of the frames of the image series together and taking the mean projection of this subset, though if the series has significant motion, the mean projection will be blurry and thus will not be a great template. A more modern solution to generating a template is to take a mean projection of a stabilized subset of frames that are well correlated to each other - this can exclude frames with significant motion artifacts. These highly correlated frames can be detected automatically, though this process can be slow. Alternatively, the user might want to manually select a number of frames where the tissue is in the desired stable position to generate the initial template

When using the `multichannel_normcorre` function, the user can specify how a template should be created by using the `channel_options.templateType` setting. The `basic` will use the 1st 200 frames to create the template. The `auto` option will use 100+ higher correlating frames to use as the template. The `manual` option will open a GUI where the user can manually specify which of the 1st 200 frames should be used to create the template.


## Original readme for the NoRMCorre repo:

# NoRMCorre: Non-Rigid Motion Correction 
This package provides a Matlab implementation of the NoRMCorre algorithm [[1]](#ref), and can be used for online piecewise rigid motion correction of 2d (planar) or 3d (volumetric) calcium imaging data. 

## Citation

If you find this package useful please cite the companion paper [[1]](#ref):

```
@article{pnevmatikakis2017normcorre,
  title={NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data},
  author={Pnevmatikakis, Eftychios A and Giovannucci, Andrea},
  journal={Journal of neuroscience methods},
  volume={291},
  pages={83--94},
  year={2017},
  publisher={Elsevier}
}
```

## Synopsis

The algorithm operates by splitting the field of view into a set of overlapping patches. For each patch and each frame a rigid translation is estimated by aligning the patch against a template using an efficient, FFT based, algorithm for subpixel registration [[2]](#reg). The estimated set of translations is further upsampled to a finer resolution to create a smooth motion field that is applied to a set of smaller overlapping patches. Extra care is taken to avoid smearing caused by interpolating overlapping patches with drastically different motion vectors. The registered frame is used to update the template in an online fashion by calculating a running/mean of past registered frames. The pipeline is summarized in the figure below.

![Alt text](pipeline.png?raw=true "piecewise rigid motion correction pipeline")

## Code details

See the function [```demo.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo.m) for an example of the code. The algorithm is implemented in the function ```normcorre.m```. If you have access to the parallel computing toolbox, then the function ```normcorre_batch.m``` can offer speed gains by enabling within mini-batch parallel processing. The user gives a dataset (either as 3D or 4D tensor loaded in RAM or memory mapped, or a pointer to a .tiff stack or .hdf5 file), and a parameters struct ```options```. Optionally, an initial template can also be given. The algorithm can also be used for motion correction of 1p micro-endoscopic data, by estimating the shifts on high pass spatially filtered version of the data. See the script [```demo_1p.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_1p.m) for an example.

The algorithm can also be ran using the MotionCorrection object. See [```demo_mc_class.m```](https://github.com/simonsfoundation/NoRMCorre/blob/master/demo_mc_class.m) for an example on how to use the object for 2p and 1p data.

The ```options``` struct can be set either manually or by using the function ```NoRMCorreSetParms.m```. Some parameters of the ```options``` struct are the following:

| Parameter name | Description |
|----------------|-------------|
| ```d1,d2,d3``` | dimensions of field of view |
| ```grid_size``` | size of non-overlapping portion of each patch the grid in each direction (x-y-z)|
| ```overlap_pre```| size of overlapping region in each direction before upsampling  |
| ```mot_uf```    | upsampling factor for smoothing and refinement of motion field |
| ```overlap_post ``` | size of overlapping region in each direction after upsampling |
| ```max_shift``` | maximum allowed shift for rigid translation | 
| ```max_dev``` | maximum deviation of each patch from estimated rigid translation |
| ```upd_template``` | update the template online after registering some frames |
| ```bin_width``` | length of bin over which the registered frames are averaged to update the template |
| ```init_batch``` | number of frames to be taken for computing initial template |
| ```iter``` | number of times to go over the dataset |
| ```output_type``` | type of output registered file |
| ```phase_flag``` | flag for using phase correlation |
| ```correct_bidir``` | check for offset due to bidirectional scanning (default: true) |

The performance of registration can be evaluated using the function ```motion_metrics.m```. The function simply computes the correlation coefficient of each (registered) frame, with the mean (registered) frame across time, the mean registered frame, and its crispness.

## Developers

[Eftychios A. Pnevmatikakis](https://github.com/epnev), Flatiron Institure, Simons Foundation

## External packages

This package includes functions from the following packages
- [Save and load a multipage tiff file](https://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-a-multiframe-tiff-image/content/loadtiff.m)
- [Savefast](https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly) for saving (and then loading) MAT files more quickly without compressing their contents. 
- [Eficient subpixel registration by cross-correlation](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) for fast alignment of an image against a template.

## Integrations 

This package will be integrated with the [Matlab code](https://www.github.com/epnev/ca_source_extraction) for source extraction and deconvolution using CNMF.

A python version of this algorithm developed from [Andrea A. Giovannuci](https://github.com/agiovann) is included as part of the [CaImAn](https://github.com/simonsfoundation/CaImAn) package that provides a complete pipeline for calcium imaging data pre-processing.

Although the two implementations give almost identical results for the same input file, there are some slight differences in the way they are called and their capabilities. These differences are highlighted [here.](https://github.com/simonsfoundation/NoRMCorre/wiki/Differences-between-Matlab-and-Python-implementations)

## More details, contact information, and citing NoRMCorre

Check the [wiki](https://github.com/simonsfoundation/NoRMCorre/wiki) for more details and some frequently asked questions. 

Please use the [gitter chat room](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) for questions and comments, and create an issue for any bugs you might encounter.

If you find this package useful please cite the following paper:

Eftychios A. Pnevmatikakis and Andrea Giovannucci, *NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data*, Journal of Neuroscience Methods, vol. 291, pp 83-94, 2017; doi: [https://doi.org/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)

## Acknowledgements

The 2p example dataset is kindly provided from Andrea Giovannucci, taken at Wang lab (Princeton University).
The 1p example dataset is kindly provided by Daniel Aharoni and Peyman Golshani (UCLA, [Miniscope project](http://miniscope.org)).

## References 

<a name="ref"></a>[1] Eftychios A. Pnevmatikakis and Andrea Giovannucci, *NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data*, Journal of Neuroscience Methods, vol. 291, pp 83-94, 2017; doi: [https://doi.org/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)

<a name="reg"></a>[2] Guizar-Sicairos, M., Thurman, S. T., & Fienup, J. R. (2008). Efficient subpixel image registration algorithms. Optics letters, 33(2), 156-158. Matlab implementation available [here](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation).
