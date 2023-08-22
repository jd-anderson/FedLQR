function cost = compute_cost(A,B,Q,R,K)

d=size(A,1);

T=1e2; %time horizon

%number of samples to compute the expectation
nsamples=10;
costs=[];
for i=1:nsamples
    x=normrnd(1,1e-6,[d,1]);
    cost_t=0;
    for t=1:T
        cost_t=cost_t + (x'*Q*x + (-K*x)'*R*(-K*x));
        x=(A-B*K)*x;
    end
    costs=[costs cost_t];
end
cost=mean(costs);
end