%%%%This is the code for ATPRK produced by Dr Qunming Wang; Email: wqm11111@126.com
%%%%Copyright belong to Qunming Wang
%%%%When using the code, please cite the fowllowing papers
%%%%Q. Wang, W. Shi, Z. Li, P. M. Atkinson. Fusion of Sentinel-2 images. Remote Sensing of Environment, 2016, 187: 241¨C252.
%%%%Q. Wang, W. Shi, P. M. Atkinson, Y. Zhao. Downscaling MODIS images with area-to-point regression kriging. Remote Sensing of Environment, 2015, 166: 191¨C204.

%%%Selected band scheme selects a 10m band with the largest CC for each 20m band
clear all;
load S2_20m;%%%20m bands in a image cube (6 bands)
load S2_10m;%%%10m bands in a image cube (4 bands)
s=2;
I_MS=S2_204_20m;
I_PAN=S2_204_10m;

w=1;
sigma=s/2;
PSFh=PSF_template(s,w,sigma);%%%Gaussian PSF
%PSFh=zeros((2*w+1)*s,(2*w+1)*s);PSFh(w*s+1:w*s+s,w*s+1:w*s+s)=1/s^2;%%%Ideal square wave PSF

%%%%correlation analysis
I_PAN_upscaled=dowmsample_cube(I_PAN,s,w,PSFh);
for i=1:6
    for j=1:4
        [RMSE0,CC0]=evaluate_relation(I_MS(:,:,i),I_PAN_upscaled(:,:,j));
        CC_matrix(i,j)=CC0;
    end
end
[II,JJ]=max(CC_matrix,[],2);

Sill_min=1;
Range_min=0.5;
L_sill=20;
L_range=20;
rate=0.1;
H=20;

tic
for i=1:6
    [xrc1,RB0,Z0]=ATPRK_PANsharpen(I_MS(:,:,i),I_PAN(:,:,JJ(i)),Sill_min,Range_min,L_sill,L_range,rate,H,w,PSFh);
    Z(:,:,i)=Z0;
end
alltime=toc

FalseColorf=Z(:,:,[3,2,1]);xf=imadjust(FalseColorf/1000,stretchlim(FalseColorf/1000),[]);figure,imshow(xf);
