
%%%     Code associated with the paper 
%%%     "On Separable Quadratic Lyapunov Functions for Convex Design of Distributed Controllers"
%%%
%%%     Authors: Luca Furieri, Yang Zheng, Antonis Papachristodoulou, Maryam
%%%     Kamgarpour
%%%
%%%     All rights reserved

%%%     To be cited as follows:
%{

@inproceedings{furieri2019separable,
  title={On separable quadratic Lyapunov functions for convex design of distributed controllers},
  author={Furieri, Luca and Zheng, Yang and Papachristodoulou, Antonis and Kamgarpour, Maryam},
  booktitle={2019 18th European Control Conference (ECC)},
  pages={42--49},
  year={2019},
  organization={IEEE}
}

%}


clear;close all


Nn = [4];
L = 16;                         %number of agents with full information

Num = length(Nn);
J   = zeros(Num,7);
degrees =zeros(1,2);
Index = Num;

%some inizialitizations
Jo1=0;
Jo2=0;
Jo1new=0;
Jo2new=0;
Jo1new2=0;
Jo2new2=0;
Jc=0;
leaders_vec=[14, 15, 3, 11, 2, 5, 9, 16, 8, 13, 7, 1, 12, 6, 10, 4];


n         = Nn(Index);
[Gp,Dist] = MeshGraph(n);  %% Plant Graph

%% Dynamics
A  = cell(n^2,n^2);   
B1 = cell(n^2,1);    

% Dynamical part
for i = 1:n^2
    B1{i} = [0;1];
    A{i,i} = [1 1; 1 2];
    for j = i+1:n^2
        if Gp(i,j) == 1
            A{i,j} = 15*exp(-norm(Dist(i,:)-Dist(j,:)))*eye(2); 
            A{j,i} = 15*exp(-norm(Dist(i,:)-Dist(j,:)))*eye(2);
        end
    end
end
B2 = B1;

[As, B1s, B2s] = NetStateModel(A,B1,B2,Gp);


%%Initial communication graphs with 0 leaders
Gc = Gp;                             % mesh topology 
Gc2=MeshGraphReduced(n);             % maximal cliques



for leaders = 0 : 1 : L
    
    
    SP = bin(kron(Gc+eye(n^2),ones(1,2)));  %% sparsity patten +eye(n^2)
    SP2= bin(kron(Gc2+eye(n^2),ones(1,2)));
    
    %% Performance
    Q  = eye(2*n^2); R = eye(1*n^2);
    
    
    
    %% Structured Optimal control P2: LMI
    
    [Ko1,Jo1,Jdiag] = StrucH2LMI(A,B1,B2,Gp,Q,R,SP);                                                % Block Diagonal approach
    [Ko1new,Jo1new,Jdiagnew,rnew] = StrucH2LMI_new(A,B1,B2,Gp,Q,R,SP);                    % Sparsity Invariance approach with T = S
    [Ko1new2,Jo1new2,Jdiagnew2,rnew2] = StrucH2LMI_new(A,B1,B2,Gp,Q,R,SP2);         %  Sparsity Invariance approach with T = T_cliques
    

    Kc = lqr(As,B2s,Q,R);
    Jc = sqrt(trace(lyap((As - B2s*Kc)',Q + Kc'*R*Kc)*(B1s*B1s')));                             %optimal centralized cost
    
  
   
    if(leaders<L)
        Gc(leaders_vec(leaders+1),:) = ones(1,n^2);
        Gc2(leaders_vec(leaders+1),:) = ones(1,n^2);
    end
    
    
    J(leaders+1,:) = [Jo1,Jo2,Jo1new,Jo2new,Jo1new2,Jo2new2,Jc];
    degrees(leaders+1,:)=[rnew,rnew2];
    
end

plots;






