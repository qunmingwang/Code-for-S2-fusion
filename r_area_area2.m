function raa=r_area_area2(h,s,xX);
Assume_L1=zeros(h+1,1);[M1,N1]=find(Assume_L1==0);
Assume_L2=zeros(s,s);[M2,N2]=find(Assume_L2==0);
for i=1:h+1
    raa(i,1)=0;
    for m=1:s^2
        for n=1:s^2
            p1=[(M1(i)-1)*s+M2(m),(N1(i)-1)*s+N2(m)];
            p2=[(M1(1)-1)*s+M2(n),(N1(1)-1)*s+N2(n)];
            raa(i,1)=raa(i,1)+myfun2(xX,norm(p1-p2));
        end
    end
end
raa=raa/s^4;