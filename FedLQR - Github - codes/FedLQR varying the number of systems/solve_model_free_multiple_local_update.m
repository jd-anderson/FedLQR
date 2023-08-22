function [cost_global_iter, Ks] = solve_model_free_multiple_local_update(A,B,R,Q,K_0,M,N,L,m,lr,r)


d=size(A{1},1);

% First we initialize all the systems with K_0

Kc=[]; %clients controller
Ks=K_0; %server controller

eta_g=2e-2;
Dels=[];
cost_initial=compute_cost(A{1},B{1},Q,R,Ks);  % Initial cost for the first system  
cost_global_iter=[cost_initial];
for n=1:N %Global iterations
     %Initialize each system with the server current controller
     for k=1:M
        Kc{k}=Ks;
     end
     for j=1:M % one local iteration for each system
         %Compute the gradient of the cost function:
         for l=1:L
             %Estimate gradient
             grad_C=estimate_gradient(A{j},B{j},Kc{j},R,Q,m,lr,r);
             eta=1e-4;
             Kc{j}=Kc{j}-eta*grad_C; % do L local iterations of GD step for each system
         end

         Dels{j}=Kc{j}-Ks; %send the local model update to the server
     end

     Del=zeros(d,d);
     for k=1:M
        Del=Del+Dels{k};
     end

     Ks=Ks+(eta_g)*(Del/M); %Aggregation
     eta_g=eta_g/1.0005;
     cost_global_iter=[cost_global_iter compute_cost(A{1},B{1},Q,R,Ks)];
end 
end