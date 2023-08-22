function cost = simulate_system(A,B,K,R,Q,l)


d=size(A,1);
x=normrnd(0,1,[d,1]);
cost=0;
for t=1:l
    x=(A-B*K)*x;
    cost=cost + (x'*Q*x + (-K*x)'*R*(-K*x));
end


end