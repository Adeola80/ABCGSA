
% clc;
% clear;
% close all;

%NOTE MODIFICATION:   1) Search Space Division 2) SCOUT BEE modification
max_it=1; 
ElitistCheck=1; 
Rpower=1;
min_flag=1; % 1: minimization, 0: maximization
F_index=1;

%% Problem Definition

CostFunction=@(x) fitnessfusion(x);        % Cost Function

nVar=10;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=-10;         % Decision Variables Lower Bound
VarMax= 10;         % Decision Variables Upper Bound

%% ABC Settings

MaxIt=100;              % Maximum Number of Iterations

nPop=100;               % Population Size (Colony Size)

nOnlooker=nPop;         % Number of Onlooker Bees

L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)

a=1;                    % Acceleration Coefficient Upper Bound

%% Initialization

% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];

% Initialize Population Array
pop=repmat(empty_bee,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Population
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% ABC Main Loop

for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    %Disruptive Selection strategy
    popcost=pop.Cost;
    if popcost(1)>0
        newpopCost=1./popcost;
    else
        newpopCost=abs(popcost);
    end
    F=zeros(nPop,1);
    [rrr ccc]=size(newpopCost);
%     see1=newpopCost
%     see2=mean([pop.Cost])
    MeanCost = newpopCost-mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        [row col]=size(newbee.Position);
N=row; 
down_min=min(min(newbee.Position,[],2));
up_max=max(max(newbee.Position,[],2));
[Fbest,Lbest,BestChart,MeanChart,wb]=GSA(F_index,N,max_it,ElitistCheck,min_flag,Rpower,down_min,up_max,newbee.Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
%         seecost=pop(i).Cost;
%         seeFbest(i)=Fbest;
        % Comparision
        if newbee.Cost<=Fbest%pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);  
            % INCLUDE NEW SCORE BEE a new equation generating the new position for scout bee 
            pop(i).Cost=CostFunction(pop(i).Position);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
    
%% Results
seebest=BestSol.Position
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
