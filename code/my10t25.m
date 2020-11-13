function  A=my10t25(Chrom,T)
global Area
b=0;
AA=Area(Chrom);
for i=1:size(T,2)
       b=b+T(i);
       a=b+1-T(i);
       c=a:b;
       A(:,c)=repmat(AA(:,i),1,T(i));
end

end
