dtf_fisher
==========

Improved Dense Trajectory features encoded by Fisher Vectors

This is the implementations of DTF+Fisher vector on UCF101 dataset.
VGG fisher vector implementation(http://www.robots.ox.ac.uk/~vgg/software/enceval_toolkit/.) and VLFfeat(http://www.vlfeat.org/) package are used. LIBSVM is also needed for one-vs-rest classification. 

For theoretical introduction about this project, please refer to posts [Action Recognition with Fisher Vectors](http://www.bo-yang.net/2014/04/30/fisher-vector-in-action-recognition) and [DTF Notes](http://www.bo-yang.net/2014/01/10/dense-trajectory-notes) on [my blog](http://bo-yang.github.io).

To run the program, please execute `master` in Matlab.

[Warning for Windows users]: Shell script [gzip_dtf_files](https://github.com/bo-yang/dtf_fisher/blob/master/gzip_dtf_files) is used to compress DTF data in `compute_fisher.m` and `subsample.m`. If you don't have Shell installed, please remove the related lines.
