function rvV=r_fine_coarse2(p_vm,W,s,xX,PSF);
Assume_L1=zeros(2*W+1,2*W+1);[M1,N1]=find(Assume_L1==0);
for i=1:(2*W+1)^2
    Tvv=zeros((2*W+1)*s,(2*W+1)*s);
    for iii=1:(2*W+1)*s %%%for every sub-pixel in the local window of current coarse pixel i
        for jjj=1:(2*W+1)*s
            p1=[(M1(i)-1-W)*s+iii,(N1(i)-1-W)*s+jjj]; %%%postion of sub-pixel in local window of coarse pixel i
            Tvv(iii,jjj)=myfun2(xX,norm(p_vm-p1)); 
        end
    end
    rvV(i,1)=sum(sum(Tvv.*PSF));
end