function P_vm=ATPK_noinform_new(s,W,S,yitaX);
[c,d]=size(S);
Simulated_part=zeros(c-2*W,d-2*W);
[M1,N1]=find(Simulated_part==0);numberM1=length(M1);
M1=M1+W;N1=N1+W;
P_vm=zeros(c*s,d*s);
for k=1:numberM1
    for i=1:s
        for j=1:s
            Local_W=S(M1(k)-W:M1(k)+W,N1(k)-W:N1(k)+W);
            co=D3_D2(yitaX(i,j,1:end-1));
            P_vm((M1(k)-1)*s+i,(N1(k)-1)*s+j)=co'*reshape(Local_W,(2*W+1)^2,1);
        end
    end
end