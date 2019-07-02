%%
% main script interaction information
clear;clc;
%%% HERE LOAD YOUR DATA, CALL IT mydata %%%%%%%%%%
load('C:\Users\dmarinaz\Dropbox\code\MI_psych_networks\ex_coll.mat');mydata=data;clear data;
%load('C:\Users\dmarinaz\Dropbox\code\MI_phys_networks\ptsd.mat');mydata=data;clear data;
%%%
[npoints, n]=size(mydata); %make sure that the variables are the 2nd dimension
p_val=0.05; %p value for surrogates
ndmax=floor(n/5); %number of variables for partial conditioning, can be changed
condtype=3; % 1 full conditioning; 2 partial conditioning; 3 triplet conditioning

%%
%%% now build the 3D matrix of II values, plus a list of red, syn,
%%% independent triplets
p_corr=p_val*6/(n*(n-1)*(n-2)); %threshold for triplets with Bonferroni
II_tot=zeros(n,n,n);
Ind_red=0;
Ind_syn=0;
Ind_ind=0;
list_red=[];
list_syn=[];
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
                list_syn=[list_syn;[i,j,k]];
            elseif Itest<0
                Ind_red=Ind_red+1;
                list_red=[list_red;[i,j,k]];
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
                list_cond=[list_red; list_syn];
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
[c,p]=corr(mydata);c=c.*(p<p_corr);
[pc,p]=partialcorr(mydata);pc=pc.*(p<p_corr);
for i=1:n;c(i,i)=0;pc(i,i)=0;end
diff_syn=zeros(n);
for i=1:n
    for j=i+1:n
        for k=j+1:n
            var_set=[i,j,k];
            if any(sum(list_syn'-var_set')==0)
                CMI_01=(CMI_binary>0);
                MI_01=(MI_binary>0);
                w=CMI_01.*MI_01;
                diff_syn=diff_syn+(CMI_01-MI_01);
                CMI_binary=CMI_binary.*w;
            end
        end
    end
end
diff_syn=triu(diff_syn);

C_plot=triu(pc)+triu(c,1)';  
MI_plot=triu(CMI_binary)+triu(MI_binary,1)';

figure;
a1=subplot(2,1,1);imagesc(C_plot,[-max(max(abs(C_plot))) max(max(abs(C_plot)))]);axis square;
title('C - lower tri: pairwise, upper tri: conditioned');
colormap(a1,brewermap([],'PRGn'));colorbar
xticks(a1,1:n);yticks(a1,1:n);
a2=subplot(2,1,2);imagesc(MI_plot);axis square;
title('MI - lower tri: pairwise, upper tri: conditioned, red: synergistic FP');
colormap(a2,brewermap([],'BuGn'));colorbar
xticks(a2,1:n);yticks(a2,1:n);
hold on
[xd, yd]=find(diff_syn);
for i=1:length(xd)
    r = rectangle('Position',[yd(i)-.5 xd(i)-.5 1 1],'EdgeColor','r','LineWidth',3);
end
set(findobj(gcf,'type','axes'),'FontSize',12)