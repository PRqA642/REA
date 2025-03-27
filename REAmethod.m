function [eT0]=et0_compute(T,w)
% Linear detrending, compute 30 years'movmean,
% et0=max(movmean)-min(movmean)
%input:target matrix T
%output:eT0 (natural variation of the input dataset).
    [a,b,t]=size(T);
    TT=reshape(T,[a*b,t]);
    T_detrend=detrend(TT);
    T_detrend_movm=movmean(T_detrend,w,2);
    et0=max(T_detrend_movm,[],2)-min(T_detrend_movm,[],2);
    eT0=reshape(et0,[a,b]);

end

function [dT,Sig_t,R_B,R_D,RR]=REA(T,m,n,et0,ts)
%% REA
% INPUT:T (target matrix),m (constant,contributes to RB),n(constant,contributes to RD),et0 (natural variation),ts(iter times)
% OUTPUT:dT (output dataset),Sig_t (quantity of uncertainty),
%%%%%%% R_B/R_D (measure of the collective model reliability with respect to the two criteria separately),RR (Weight matrix)
    [a,b,c,k]=size(T);
    B=zeros(a,b,c,k,'single');

    D=zeros(a,b,c,k,'single');
    dT0=mean(T,4);
    dT1=dT0;
    for j=1:ts
        S1=0;
        S2=0;
        S3=0;
        RR=zeros(a,b,c,k,'single');
        for i=1:k
            B(:,:,:,i)=(T(:,:,:,i)-dT0);
            D(:,:,:,i)=(T(:,:,:,i)-dT1);
            Btemp=B(:,:,:,i);
            Dtemp=D(:,:,:,i);
            tempb=(abs(Btemp)<1);
            Btemp(tempb)=1;
            tempd=abs(Dtemp)<1;
            Dtemp(tempd)=1;
            B(:,:,:,i)=Btemp;
            D(:,:,:,i)=Dtemp;
            RR(:,:,:,i)=(((et0(:,:,i)./abs(B(:,:,:,i))).^m).*((et0(:,:,i)./abs(D(:,:,:,i))).^n)).^(1/(m*n));
            RR(T==0)=0;
            RT=RR(:,:,:,i).*T(:,:,:,i);
            S1=S1+RT;
            S2=S2+RR(:,:,:,i);
            S3=S3+RR(:,:,:,i).^2;
        end
        dT1=S1./S2;
    end
    dT=dT1;
    S4=0;
    for i=1:k    
        S4=S4+RR(:,i).*((T(:,:,:,i)-dT).^2);
    end
    Sig_t=(single(S4)./single(S2)).^(1/2);
    SB=sum(B,4);
    R_B=single(SB)./k;
    SD=sum(D,4);
    R_D=single(SD)./k;
end
