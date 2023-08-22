function [A,B] = sysgen(A_0,B_0,V,U,M,eps_1,eps_2,rand_n)

A={};
B={};
for i=1:M(end)
   if i==1
       A{i}=A_0;
       B{i}=B_0;
   else
       A{i}=A_0+eps_1*rand_n(i)*V;
       B{i}=B_0+eps_2*rand_n(i)*U;
   end
end 

end