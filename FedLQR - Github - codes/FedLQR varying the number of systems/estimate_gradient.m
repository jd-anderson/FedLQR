function grad_C=estimate_gradient(A,B,K,R,Q,m,l,r)


d=size(A,1);%dimension d
k=size(B,2);%dimension k
sigu = 0.5; %standard deviation of the matrices Uj

Ks=[];
cost_emp=[];
Us=[];
for j=1:m
    rand_matrix=normrnd(0,sigu,[d,d]);
    Us{j}=(r*(rand_matrix))/norm(rand_matrix,"fro");
    %Sample policy
    Ks{j}=K+Us{j};
    %Compute empirical cost through simulations
    cost_emp{j}=simulate_system(A,B,Ks{j},R,Q,l);
end

sum_cost=zeros(d,d);
for j=1:m
    sum_cost=sum_cost + cost_emp{j}*((d*k)/(r^2))*Us{j};
end

grad_C=(1/m)*sum_cost;

end