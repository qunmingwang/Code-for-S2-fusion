function rh=semivariogram(J,h);
[a,b]=size(J);
N1=0;r1=0;
for i=h+1:a
    for j=1:b
        r1=r1+(J(i,j)-J(i-h,j))^2;
        N1=N1+1;
    end
end
N2=0;r2=0;
for i=1:a-h
    for j=1:b
        r2=r2+(J(i,j)-J(i+h,j))^2;
        N2=N2+1;
    end
end
N3=0;r3=0;
for i=1:a
    for j=h+1:b
        r3=r3+(J(i,j)-J(i,j-h))^2;
        N3=N3+1;
    end
end
N4=0;r4=0;
for i=1:a
    for j=1:b-h
        r4=r4+(J(i,j)-J(i,j+h))^2;
        N4=N4+1;
    end
end
r=r1+r2+r3+r4;
N=N1+N2+N3+N4;
rh=r/(2*N);