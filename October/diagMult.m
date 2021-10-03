function product = diagMult(mtr, diag)
    [m, n]=size(mtr);
    product=zeros([m,n]);
    for i=1:n
        product(:,i)=diag(i)*mtr(:,i);
    end
end
