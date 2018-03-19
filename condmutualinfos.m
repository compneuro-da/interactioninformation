function [h, h0] = condmutualinfos(vec1,vec2,condvec,th)

nperm=1000; % number of permutations

h0 =condmutualinfo(vec1,vec2,condvec);
n=length(vec1);
ncont=0;
if h0==0
    h=0;
else
    for j=1:nperm
        vec1=vec1(randperm(n));
        h=condmutualinfo(vec1,vec2,condvec);
        if h>h0
            ncont=ncont+1;
        end
    end
    if ncont >th*nperm
        h=0;
    else
        h=h0;
    end
end