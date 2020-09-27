function TVV=T_coarse_coarse2(W,s,xX,PSF);
Assume_L1=zeros(2*W+1,2*W+1);[M1,N1]=find(Assume_L1==0);
TVV=zeros((2*W+1)^2,(2*W+1)^2);
for i=1:(2*W+1)^2
    for j=1:(2*W+1)^2 %%%for every coarse pixel i and j
        TvV=zeros((2*W+1)*s,(2*W+1)*s);
        for ii=1:(2*W+1)*s
            for jj=1:(2*W+1)*s %%%for every sub-pixel in the local window of coarse pixel j 
                
                Tvv=zeros((2*W+1)*s,(2*W+1)*s);
                for iii=1:(2*W+1)*s %%%for every sub-pixel in the local window of current coarse pixel i
                    for jjj=1:(2*W+1)*s
                        p1=[(M1(i)-1-W)*s+iii,(N1(i)-1-W)*s+jjj]; %%%postion of sub-pixel in local window of coarse pixel i
                        p2=[(M1(j)-1-W)*s+ii,(N1(j)-1-W)*s+jj]; %%%postion of sub-pixel in local window of coarse pixel j 
                        Tvv(iii,jjj)=myfun2(xX,norm(p1-p2));               
                    end
                end
                TvV(ii,jj)=sum(sum(Tvv.*PSF));
                
            end
        end
        
        TVV(i,j)=sum(sum(TvV.*PSF));
    end
end