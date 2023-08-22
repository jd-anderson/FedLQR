% ==========================================================================
% FedLQR: Federated-based policy learning for the LQR problem
% Han Wang, Leonardo F. Toso, Aritra Mitra,  James Anderson
% ==========================================================================

clc;clear all; close all

%% System Data:

% System dimensions

d=3; % nx - number of states
k=3; % nu - number of inputs

% Nominal system: 

A_0=[1.2    0.5   0.4
     0.01   0.75  0.3
     0.1    0.02  1.5];

B_0=eye(k);

V=eye(d);U=eye(d);

Q=2*eye(d);

R=0.5*eye(d);

K_0=1.62*eye(d);

%Number of clients:

M_vec=[1 10 50];


%Dissimilarity:

eps=1e-2;
eps_1=eps;
eps_2=eps;


%Gradient estimation - Zeroth order estimation parameters

r=0.1;
m=5;%number of trajectories
l=15; %rollout length
 
% Generating the system matrices

[A,B] = sysgen(A_0,B_0,V,U,M_vec,eps_1,eps_2); %This function generates similar systems
   
cost_initial=compute_cost(A{1},B{1},Q,R,K_0);

%Compute the optimal controller for the system 1

K_opt_1=dlqr(A{1},B{1},Q,R);  %Optimal control gain related to the first system
cost_opt_1_init=compute_cost(A{1},B{1},Q,R,K_opt_1); %Compute the optimal cost for this optimal control gain,


N=0.5e4; %number of global iterations
L=1;%number of local iterations
nr=20;%number of realizations
costs=[];
[A,B] = sysgen(A_0,B_0,V,U,M_vec,eps_1,eps_2); %This function generates similar systems
K_opt_1=dlqr(A{1},B{1},Q,R);
cost_opt_1=compute_cost(A{1},B{1},Q,R,K_opt_1);

norm_1=[];
norm_2=[];
for p=1:M_vec(end)
    for s=1:M_vec(end)
        norm_1=[norm_1 norm(A{p}-A{s})];
        norm_2=[norm_2 norm(B{p}-B{s})];
    end
end


dis1=["eps1=",max(norm_1)];
dis2=["eps2=",max(norm_2)];

disp(dis1);
disp(dis2);

for j=1:nr
    costs_M=zeros(length(M_vec),N+1);
    for i=1:length(M_vec)
         M=M_vec(i);
         A_matrices=[];
         B_matrices=[];
         for k=1:M
             A_matrices{k}=A{k};
             B_matrices{k}=B{k};
         end
          % Model-free - Federated Learning
          [costs_M(i,:),K] = solve_model_free_multiple_local_update(A_matrices,B_matrices,R,Q,K_0,M,N,L,m,l,r);
          costs_M(i,:)=(costs_M(i,:) - cost_opt_1)/cost_opt_1;
    end
    costs{j}=costs_M;
end