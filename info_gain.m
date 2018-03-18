function [y, ind]=info_gain(drive,X,nvar,ndmax)
indx=setdiff(1:nvar,drive); %eliminate the candidate driver from the set
t=X{drive};
% t=X(:,:,drive);
Zt=[];
for nd=1:ndmax
    n1=length(indx);
    z=zeros(n1,1);
    for k=1:n1
        Zd=mergemultivariables(Zt, X{indx(k)});
        z(k)= mutualinfo(Zd,t); %compute MI
    end
    [y(1,nd), id]=max(z); %greedy algorithm, find the max contribution, store it and remove it from the set of candidates
    Zt=mergemultivariables(Zt, X{indx(id)});
    ind(1,nd)=indx(id);
    indx=setdiff(indx,indx(id));
end