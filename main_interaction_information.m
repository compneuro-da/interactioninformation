%%
% main script interaction information
clear;clc;
%%% HERE LOAD YOUR DATA, CALL IT mydata %%%%%%%%%%
load('C:\Users\dmarinaz\Dropbox\code\MI_phys_networks\ptsd.mat');mydata=data;clear data;
%%%
[npoints, n]=size(mydata); %make sure that the variables are the 2nd dimension
p_val=0.05; %p value for surrogates
ndmax=floor(n/10); %number of variables for partial conditioning
condtype=3; % 1 full conditioning; 2 partial conditioning; 3 triplet conditioning

%%
%%% now build the 3D matrix of II values, plus a list of red, syn,
%%% independent triplets
p_corr=p_val*6/(n*(n-1)*(n-2)); %threshold for triplets with Bonferroni
II_tot=zeros(n,n,n);
Ind_red=0;
Ind_syn=0;
Ind_ind=0;
for i=1:n
    for j=i+1:n
        for k=j+1:n
            [Itest, II]=interaction_inf(mydata(:,i),mydata(:,j),mydata(:,k),p_corr);
            II_tot(i,j,k)=II;
            II_tot(i,k,j)=II;
            II_tot(j,i,k)=II;
            II_tot(j,k,i)=II;
            II_tot(k,i,j)=II;
            II_tot(k,j,i)=II;
            if Itest>0
                Ind_syn=Ind_syn+1;
                list_syn(Ind_syn,:)=[i,j,k];
            elseif Itest<0
                Ind_red=Ind_red+1;
                list_red(Ind_red,:)=[i,j,k];
            else
                Ind_ind=Ind_ind+1;
                list_ind(Ind_ind,:)=[i,j,k];
            end
        end
    end
end

%%
p_corr=p_val/(n*(n-1)*0.5); % threshold with Bonferroni
MI_binary=zeros(n);
MI=zeros(n);
CMI_binary=MI_binary;
CMI=MI;
NET_II=MI;
for i=1:n
    for j=i+1:n
        [MI_binary(i,j), MI(i,j)]=mutualinfos(mydata(:,i),mydata(:,j),p_corr); %mutual info with threshold
        MI_binary(j,i)=MI_binary(i,j);
        MI(j,i)=MI(i,j);
        switch condtype
            case 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % here you compute MI conditioned to the rest of the system. Can be tricky with many variables/fewer points
                condind=setdiff(1:n,[i,j]);
                condvec=mydata(:,condind(1));
                for icond=2:length(condind)
                    condvec=mergemultivariables(condvec,mydata(:,condind(icond)));
                end
                [CMI_binary(i,j), CMI(i,j)]=condmutualinfos(mydata(:,i),mydata(:,j),condvec,p_corr); %mutual info with threshold
                CMI_binary(j,i)=CMI_binary(i,j);
                CMI(j,i)=CMI(i,j);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 2
                %%% here you condition to ndmax most informative variables for each
                %%% driver
                [INFO, ind_PC]=init_partial_conditioning(mydata,ndmax);
                A=ind_PC(j,:);
                condind = A(~ismembc(A(:), i));
                %condind = condind(1:ndmax);
                condvec=mydata(:,condind(1));
                for icond=2:length(condind)
                    condvec=mergemultivariables(condvec,mydata(:,condind(icond)));
                end
                [CMI_binary(i,j), CMI(i,j)]=condmutualinfos(mydata(:,i),mydata(:,j),condvec,p_corr); %mutual info with threshold
                CMI_binary(j,i)=CMI_binary(i,j);
                CMI(j,i)=CMI(i,j);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 3
                %%% here you condition only on third members of synergetic or redundant triplets
                list_cond=[list_red; list_red];
                row = find(any(list_cond == i, 2) & any(list_cond == j, 2));
                NET_II(i,j)=length(row); %build a graph in which the links are the number of multiplets in which the two nodes are both present
                NET_II(j,i)=NET_II(j,i);
                if ~isempty(row)
                    condvec=zeros(npoints,1);
                    for icond=1:length(row)
                        condind=setdiff(list_cond(row(icond),:),[i,j]);
                        condvec=mergemultivariables(condvec,mydata(:,condind));
                    end
                    [CMI_binary(i,j), CMI(i,j)]=condmutualinfos(mydata(:,i),mydata(:,j),condvec,p_corr); %mutual info with threshold
                    CMI_binary(j,i)=CMI_binary(i,j);
                    CMI(j,i)=CMI(i,j);
                else
                    CMI_binary(i,j)=MI_binary(i,j);
                    CMI_binary(j,i)=MI_binary(j,i);
                    CMI(i,j)=MI(i,j);
                    CMI(j,i)=MI(j,i);
                end
        end
    end
end