# PyramidTextureFiltering
Pyramid Texture Filtering C++ implementation
https://arxiv.org/pdf/2305.06525

**Technique and Importance of the paper**
The technique proposed in the paper follows the simple insight that the highest level in the Gaussian pyramid construction of an image erases textures while preserving the important structures. It involves upsampling the coarsest Gaussian pyramid using the next finer Gaussian level as a guide, and adding it to the corresponding Laplacian pyramid level. The output is again smoothed and is used in the upsampling in the next iteration.

The difficult aspect of texture smoothing lies in distinguishing texture from structure while smoothing, and therefore edge preserving smoothing techniques are not viable since they do not distinguish between the two. Structure aware upsampling enables this technique to produce a full resolution smoothed image, while still maintaing sharp features and avoiding visual artifacts. Another aspect of the proposal of this paper worth mentioning is that it only uses filtering techniques and does not rely on machine learning to extract image structures. 

**Implemented code**
Bilateral Filtering, Joint Bilateral Filtering, Gaussian Smoothing, Bilinear Interpolation and up/down sampling using the former two are implemented as methods. Gaussian Pyramids, Laplacian pyramids are constructed using these, to finally perform the PSU iterations.
The authors' suggestions for the ranges of hyperparameters are followed.

![Original Image](https://github.com/singh-todai/PyramidTextureFiltering/blob/main/image3.png)
            
                            |
                            |
                            V
![Smoothed Image](https://github.com/singh-todai/PyramidTextureFiltering/blob/main/output3.png)

It can be seen that the mosaic effect is nicely removed.
