function h = mutualinfos(vec1,vec2,th)
h0 = mutualinfo(vec1,vec2);
n=length(vec1);
ncont=0;
for j=1:1000
vec1=vec1(randperm(n));
h=mutualinfo(vec1,vec2); 
if h>h0
    ncont=ncont+1;
end
end
if ncont >th*1000
    h=0;
else
    h=1;
end



