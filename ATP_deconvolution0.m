%%%%%x_area is the coarse semivariogram derived from coarse proportion image
%%%%%H is the length of lags
function x_best=ATP_deconvolution0(H,s,x_area,Sill_min,Range_min,L_sill,L_range,rate);
Fa0=myfun2(x_area,[1:1:s*H]);
Fa0_vector=Fa0(s:s:end);Fa0_vector=Fa0_vector';
%[x1,resnorm]=lsqcurvefit(@myfun2,x0,xdata,rh(i,:));
Dif_min=10^6;
for i=1:L_sill%%%%sill
    for j=1:L_range%%%%%range
        xp=[(Sill_min+rate*i)*x_area(1),(Range_min+rate*j)*x_area(2)];
        raa0=r_area_area2(H,s,xp);
        raa=raa0(2:H+1,1)-raa0(1,1);
        Dif=norm(raa-Fa0_vector);
        if Dif<=Dif_min
            x_best=xp;Dif_min=Dif;
        end
    end
end