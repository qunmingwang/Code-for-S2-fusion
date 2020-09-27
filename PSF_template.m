%%%2*w+1 is the width of the PSF; b is the kernel size;s is the DS zoom factor
function H=PSF_template(s,w,b);
for i=1:(2*w+1)*s
    for j=1:(2*w+1)*s
        Dis2=(norm([i-0.5,j-0.5]-[(2*w+1)*s/2,(2*w+1)*s/2]))^2;
        H0(i,j)=exp(-Dis2/(2*b^2));
    end
end
Hsum=sum(sum(H0));
H=H0/Hsum;