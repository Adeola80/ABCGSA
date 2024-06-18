% Gravitational Search Algorithm.
function [Fbest,Lbest,BestChart,MeanChart,wb]=GSA(F_index,N,max_it,ElitistCheck,min_flag,Rpower,down_min,up_max,feat)

%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of the test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm in eq.8.
%Rpower: Power of R in eq.7.
 Rnorm=2; 
%get allowable range and dimension of the test function.
[low,up,dim]=objective_functions_range(F_index,down_min,up_max);
dd=randi([5 5]);

%random initialization for agents.
seesee=feat;
X=feat(:,1:dd);%initialization(feat,dim,N,up,low)

%create the best so far chart and average fitnesses chart.
BestChart=[];MeanChart=[];

V=zeros(N,dim);

for iteration=1:max_it
%     iteration
    
    %Checking allowable range. 
    X=space_bound(X,up,low); 

    %Evaluation of agents. 
    fitness=evaluateF(X,F_index); 
    
    if min_flag==1
    [best best_X]=min(fitness); %minimization.
    else
    [best best_X]=max(fitness); %maximization.
    end        
    
    if iteration==1
       Fbest=best;Lbest=X(best_X,:);
    end
    if min_flag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(best_X,:);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(best_X,:);
      end
    end
      
BestChart=[BestChart Fbest];
MeanChart=[MeanChart mean(fitness)];

%Calculation of M. eq.14-20
[M]=massCalculation(fitness,min_flag); 

%Calculation of Gravitational constant. eq.13.
G=Gconstant(iteration,max_it); 

%Calculation of accelaration in gravitational field. eq.7-10,21.
a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);

%Agent movement. eq.11-12
[X,V]=move(X,a,V);
Fbest=Fbest/2;
Lbest=Lbest/2;

if Fbest==0
    Fbest=(unifrnd(min(min(feat,[],2)),min(min(feat,[],2))+((min(min(feat,[],2)))/4),1,1))/2;
else
    Fbest=min(min(feat,[],2))/2;
end %iteration

if Lbest(1)>=0
    Lbest=(unifrnd(max(max(feat,[],2)),max(max(feat,[],2))+((max(max(feat,[],2)))/4),1,1))/3;
else
    Lbest=max(max(feat,[],2))-max(max(feat,[],2))/3;
end    
    wb=randperm(5);
end %iteration
