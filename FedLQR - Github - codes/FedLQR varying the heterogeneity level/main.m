% ==========================================================================
% FedLQR: Federated-based policy learning for the LQR problem
% Han Wang, Leonardo F. Toso, Aritra Mitra,  James Anderson
% ==========================================================================

clc;clear all; close all

%% System Data:


% System dimensions

d=3;
k=3;

% Nominal system: 

A_0=[1.2    0.5   0.4
     0.01   0.75  0.3
     0.1    0.02  1.5];

B_0=eye(k);


V=diag([3.5 1 0.1]);U=diag([1.5 0.1 1]);

Q=2*eye(d);

R=0.5*eye(d);

K_0=1.62*eye(d);


%Number of clients:

M_vec=[10];


%Dissimilarity:

eps=1e-2;
eps_1=eps;
eps_2=eps;


%Gradient estimation - Zeroth order estimation

r=0.1;
m=5;%number of trajectories
l=15; %rollout length


N=0.5e4; %number of global iterations
L=1;%number of local iterations
nr=20;%number of realizations
rand_n=rand(1,M_vec(end));

%Zero-heterogeneity

[A_1,B_1] = sysgen(A_0,B_0,V,U,M_vec,0,0,rand_n); %This function generates similar systems
K_opt_1=dlqr(A_1{1},B_1{1},Q,R);
cost_opt_1=compute_cost(A_1{1},B_1{1},Q,R,K_opt_1);

norm_1=[];
norm_2=[];
for p=1:M_vec(end)
    for s=1:M_vec(end)
        norm_1=[norm_1 norm(A_1{p}-A_1{s})];
        norm_2=[norm_2 norm(B_1{p}-B_1{s})];
    end
end


dis1=["eps1=",max(norm_1)];
dis2=["eps2=",max(norm_2)];

disp(dis1);
disp(dis2);
costs_1=[];

for j=1:nr
    costs_M=zeros(length(M_vec),N+1);
    for i=1:length(M_vec)
         M=M_vec(i);
         A_matrices=[];
         B_matrices=[];
         for k=1:M
             A_matrices{k}=A_1{k};
             B_matrices{k}=B_1{k};
         end
          % Model-free - Federated Learning
          [costs_M(i,:),K] = solve_model_free_multiple_local_update(A_matrices,B_matrices,R,Q,K_0,M,N,L,m,l,r);
          costs_M(i,:)=(costs_M(i,:) - cost_opt_1)/cost_opt_1;
    end
    costs_1{j}=costs_M;
end



%eps_1 = 1e-2, eps_2=1e-2

[A_2,B_2] = sysgen(A_0,B_0,V,U,M_vec,1e-1,1e-1,rand_n); %This function generates similar systems
K_opt_1=dlqr(A_2{1},B_2{1},Q,R);
cost_opt_1=compute_cost(A_2{1},B_2{1},Q,R,K_opt_1);

norm_1=[];
norm_2=[];
for p=1:M_vec(end)
    for s=1:M_vec(end)
        norm_1=[norm_1 norm(A_2{p}-A_2{s})];
        norm_2=[norm_2 norm(B_2{p}-B_2{s})];
    end
end


dis3=["eps1=",max(norm_1)];
dis4=["eps2=",max(norm_2)];

disp(dis3);
disp(dis4);
costs_2=[];

for j=1:nr
    costs_M=zeros(length(M_vec),N+1);
    for i=1:length(M_vec)
         M=M_vec(i);
         A_matrices=[];
         B_matrices=[];
         for k=1:M
             A_matrices{k}=A_2{k};
             B_matrices{k}=B_2{k};
         end
          % Model-free - Federated Learning
          [costs_M(i,:),K] = solve_model_free_multiple_local_update(A_matrices,B_matrices,R,Q,K_0,M,N,L,m,l,r);
          costs_M(i,:)=(costs_M(i,:) - cost_opt_1)/cost_opt_1;
    end
    costs_2{j}=costs_M;
end




%eps_1 = 1e-2, eps_2=1e-2

[A_3,B_3] = sysgen(A_0,B_0,V,U,M_vec,5e-1,5e-1,rand_n); %This function generates similar systems
K_opt_1=dlqr(A_3{1},B_3{1},Q,R);
cost_opt_1=compute_cost(A_3{1},B_3{1},Q,R,K_opt_1);

norm_1=[];
norm_2=[];
for p=1:M_vec(end)
    for s=1:M_vec(end)
        norm_1=[norm_1 norm(A_3{p}-A_3{s})];
        norm_2=[norm_2 norm(B_3{p}-B_3{s})];
    end
end


dis5=["eps1=",max(norm_1)];
dis6=["eps2=",max(norm_2)];

disp(dis5);
disp(dis6);
costs_3=[];

for j=1:nr
    costs_M=zeros(length(M_vec),N+1);
    for i=1:length(M_vec)
         M=M_vec(i);
         A_matrices=[];
         B_matrices=[];
         for k=1:M
             A_matrices{k}=A_3{k};
             B_matrices{k}=B_3{k};
         end
          % Model-free - Federated Learning
          [costs_M(i,:),K] = solve_model_free_multiple_local_update(A_matrices,B_matrices,R,Q,K_0,M,N,L,m,l,r);
          costs_M(i,:)=(costs_M(i,:) - cost_opt_1)/cost_opt_1;
    end
    costs_3{j}=costs_M;
end


gap_hey_1=zeros(nr,N+1);
gap_hey_2=zeros(nr,N+1);
gap_hey_3=zeros(nr,N+1);


for i=1:nr
    gap_het_1(i,:)=costs_1{i};
end 

for i=1:nr
    gap_het_2(i,:)=costs_2{i};
end 

for i=1:nr
    gap_het_3(i,:)=costs_3{i};
end 

save('gap_het_1.mat','gap_het_1');
save('gap_het_2.mat','gap_het_2');
save('gap_het_3.mat','gap_het_3');


