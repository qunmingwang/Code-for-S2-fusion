function yita=ATPK_noinform_yita_new(s,W,xX,PSF);
TVV=T_coarse_coarse2(W,s,xX,PSF);
for i=1:s
    for j=1:s
        cordinate_vm=[W*s+i,W*s+j];
        %rvV=r_fine_coarse2(cordinate_vm,W,s,xX);
        rvV=r_fine_coarse2(cordinate_vm,W,s,xX,PSF);
        Matrix=[TVV,ones((2*W+1)^2,1);ones(1,(2*W+1)^2),0];
        Vector=[rvV;1];
        yita(i,j,:)=inv(Matrix)*Vector;
    end
end