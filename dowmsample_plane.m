function S=dowmsample_plane(plane,s,w,PSF);
plane=Extend_plane(plane,w*s);
[sizec,sized]=size(plane);

S=zeros(sizec/s,sized/s);
for i=w*s+1:s:sizec-w*s
    for j=w*s+1:s:sized-w*s
        m=(i+s-1)/s; n=(j+s-1)/s;
        LW=plane(i-w*s:i+w*s+s-1,j-w*s:j+w*s+s-1);%%%Local Window with (2*w+1)*s sub-pixels
        S(m,n)=sum(sum(LW.*PSF));
    end
end
S=S(w+1:end-w,w+1:end-w);