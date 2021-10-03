function output= new_SBL_algo (sigma, phi, y, epsilon)
%tic
    [m,n]=size(phi);
    phiT=phi.';
    phiTy=phiT*y;

    reqsigma=1/(sigma*sigma);
    mu=zeros([n,1]);
    mymat=ones([n,1]);
    k=100;
    s=0;
    %toc
    
    while k>epsilon 
        oldmu = mu;
        
        tau=diag(mymat);
        
        %tic
        %phitau=phi*tau;
        phitau=diagMult(phi, mymat);
        phitauT=phitau.';
        %toc
        
        %tic
        g = sigma*sigma*eye(m)+phitau*phiT;
        %toc
        
        
        
        % inv(A)*b is A\b
        %tic
        SIGMA = tau - phitauT*(g\(phitau));
        %toc
        
        
        %tic
        mu=reqsigma*(SIGMA*phiTy); 
        %toc
        
        
        for i=1:n
            mymat(i)=mu(i)*mu(i)+SIGMA(i,i);
        end
        
        s=s+1
        %norm(y-phi*mu);
        k=norm(oldmu-mu)
        
        if s>500
            break
        end
    end
    output=mu;
end
