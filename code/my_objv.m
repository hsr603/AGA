function [F,U,ObjV]=my_objv(D)

global L p F_limit U_limit
ObjV=zeros(1,6);

[F,U]=ST(D); 
F=F./D.A'*1e-3;  
ObjV(1)=sum((abs(F)>=F_limit));
f=find((abs(F))>F_limit);
if f>0;ObjV(2)=sum(abs(F(f))/F_limit);end
ObjV(3)=sum(sum((abs(U)>=U_limit)));
f=find((abs(U))>U_limit);
if f>0;ObjV(4)=sum(abs(U(f))/U_limit);end
ObjV(5)=sum(p*D.A.*L');
k=abs(F)>F_limit;
k1=sum(k);  
k=abs(U)>U_limit;
k2=sum(sum(k)); 
k1=~k1;
k2=~k2;
ObjV(6)=(ObjV(2)+ObjV(4))*10000+ObjV(:,5)*(1*k1*k2);
end
