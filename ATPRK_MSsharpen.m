function [xrc1,RB,Z]=ATPRK_MSsharpen(Coarse,PAN,Sill_min,Range_min,L_sill,L_range,rate,H,w,PSF);
[a1,b1]=size(Coarse);
[a2,b2,c2]=size(PAN);
s=a2/a1;

%%%%%linear regression modeling
PAN_upscaled=dowmsample_cube(PAN,s,w,PSF);
GB0=D3_D2(PAN_upscaled);
GB1=reshape(Coarse,1,a1*b1);
xrc1=lsqlin([ones(a1*b1,1),GB0'],GB1,[],[]) ;
GBF=D3_D2(PAN);
Ff1=[ones(a2*b2,1),GBF']*xrc1;Z_R=reshape(Ff1,a1*s,b1*s);

%%%%%residual calculation
Z_R_upscaled=dowmsample_plane(Z_R,s,w,PSF);
RB=Coarse-Z_R_upscaled;

%%%%%ATPK for residuals, Deconvolution is achieved by a trail-and-error procedure
W=w;
RB_extend1=[repmat(RB(:,1),[1,W]),RB,repmat(RB(:,end),[1,W])];%%%extend columns
RB_extend=[repmat(RB_extend1(1,:),[W,1]);RB_extend1;repmat(RB_extend1(end,:),[W,1])];%%%extend rows

x0=[0.1,1];%%%%%x0 is the initial value for fitting
for h=1:H
    rh(h)=semivariogram(RB,h);
end
[xa1,resnorm]=lsqcurvefit(@myfun2,x0,s:s:s*H,rh);Fa1=myfun2(xa1,1:1:s*H);
xp_best=ATP_deconvolution0(H,s,xa1,Sill_min,Range_min,L_sill,L_range,rate);Fp=myfun2(xp_best,[1:1:s*H]);
raa0=r_area_area2(H,s,xp_best);raa=raa0(2:H+1,1)-raa0(1,1);[xa2,resnorm]=lsqcurvefit(@myfun2,x0,s:s:s*H,raa');Fa2=myfun2(xa2,[1:1:s*H]);
xp_best_matrix=xp_best;

yita1=ATPK_noinform_yita_new(s,W,xp_best,PSF);
P_vm=ATPK_noinform_new(s,W,RB_extend,yita1);
Z_ATPK=P_vm(W*s+1:end-W*s,W*s+1:end-W*s);
Z=Z_R+Z_ATPK;
