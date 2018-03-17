clear;clc;
%%% HERE LOAD YOUR DATA, CALL IT mydata %%%%%%%%%%
[npoints, n]=size(mydata); %make sure that the variables are the 2nd dimension
th=0.05/(n*(n-1)*0.5); % threshold with Bonferroni
for i=1:n
    for j=i+1:n
        muti(i,j)=mutualinfos(mydata(:,i),mydata(:,j),th); %mutual info with threshold
        muti(j,i)=muti(i,j);
    end
end

th=0.05*6/(n*(n-1)*(n-2)); %threshold for triplets with Bonferroni

h=0;
for i=1:n
    for j=i+1:n
        for k=j+1:n
            h=h+1;
            tripl(h)=interaction_inf(mydata(:,i),mydata(:,j),mydata(:,k),th);
        end
    end
end


