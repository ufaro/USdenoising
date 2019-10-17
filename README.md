# USdenoising

USdenoising by Younes Farouj. v.April 2017


Ultrasound spatiotemporal despeckling. This Matlab code follows the approaches presented in the two papers; 

[1] Y. Farouj, L. Navarro, J.-M. Freyermuth, M. Clausel, and P. Delachartre. Ultrasound Spatio- temporal Despeckling via Kronecker Wavelet-Fisz Thresholding. Signal, Image and Video Processing, 2018, vol. 12, no 6, p. 1125-1132

and 

[2] Y. Farouj, J.-M. Freyermuth, L. Navarro, M. Clausel, and P. Delachartre. Hyperbolic Wavelet-Fisz denoising for a model arising in Ultrasound Imaging. IEEE Transactions on Computational Imaging, vol. 3, no. 1, pp. 1–10, 2017. 

The code includes two version. A data-driven one and parametric one. 



INPUT

  
% %    Finds the data-driven wavelet-Fisz estimate using some given (BseOnd) wavelet.
% %    
% %    takes:
% %    x           - data: An ultrasound imaging sequence  
% %    h           - variance function (for the parametric version)
% %    AV          - a pre-estimate of the distribution of pixel values. Can be an
% %                  average overt time of the dynamic sequence. (for the data-driven version)
% %    zero.levels - number of finest-scale levels at which detail coefficients are set to       zero 
% %    thr         - type of threshold: 1 for universal, 2 set by user
% %    theta       - multiplier for the (universal) threshold (use 1 if not sure)
% %    zero.levels - number of finest-scale levels at which detail coefficients are set to zero
% %    BseOnd      - wavelet type. 

OUTPUT

Reconstructed: filtered sequence

Parameters for the algorithm (can be changed by user)

If any problem, please contact Younes Farouj via younes.farouj@epfl.ch.

