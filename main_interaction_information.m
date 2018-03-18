%%
% main script interaction information
clear;clc;
%%% HERE LOAD YOUR DATA, CALL IT mydata %%%%%%%%%%
load('C:\Users\dmarinaz\Dropbox\code\MI_phys_networks\santos.mat');mydata=data;clear data;
[npoints, n]=size(mydata); %make sure that the variables are the 2nd dimension
th=0.05/(n*(n-1)*0.5); % threshold with Bonferroni

%%
MI_binary=zeros(n);
MI=zeros(n);
CMI_binary=zeros(n);
CMI=zeros(n);
for i=1:n
    for j=i+1:n
        [MI_binary(i,j), MI(i,j)]=mutualinfos(mydata(:,i),mydata(:,j),th); %mutual info with threshold
        MI_binary(j,i)=MI_binary(i,j);
        MI(j,i)=MI(i,j);
        % here you compute MI conditioned to the rest of the system. Can be tricky with many variables/fewer points 
        condind=setdiff(1:n,[i,j]);
        condvec=zeros(npoints,1);
        for icond=1:n-2
            condvec=mergemultivariables(condvec,mydata(:,condind(icond)));
        end
        [CMI_binary(i,j), CMI(i,j)]=condmutualinfos(mydata(:,i),mydata(:,j),condvec,th); %mutual info with threshold
        CMI_binary(j,i)=CMI_binary(i,j);
        CMI(j,i)=CMI(i,j);
    end
end

th=0.05*6/(n*(n-1)*(n-2)); %threshold for triplets with Bonferroni
%%
%%% now build the 3D matrix of II values, plus a list of red, syn,
%%% independent triplets
II_tot=zeros(n,n,n);
Ind_red=0;
Ind_syn=0;
Ind_ind=0;
for i=1:n
    for j=i+1:n
        for k=j+1:n
            [Itest, II]=interaction_inf(mydata(:,i),mydata(:,j),mydata(:,k),th);
            II_tot(i,j,k)=II;
            II_tot(i,k,j)=II;
            II_tot(j,i,k)=II;
            II_tot(j,k,i)=II;
            II_tot(k,i,j)=II;
            II_tot(k,j,i)=II;
            if Itest>0
                Ind_syn=Ind_syn+1;
                list_syn{Ind_syn}=mat2str([i,j,k]);
            elseif Itest<0
                Ind_red=Ind_red+1;
                list_red{Ind_red}=mat2str([i,j,k]);
            else
                Ind_ind=Ind_ind+1;
                list_ind{Ind_ind}=mat2str([i,j,k]);
            end
        end
    end
end