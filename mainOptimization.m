%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      O P T I M I Z A T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         This script optimizes the model and saves the results.
%         It solves the model using various optimization methods.
%
%         The user must choose between scenarios 1, 2, and 3 (models) from the paper
%         or introduce a new one, for example, by modifying one of the
%         existing ones.
%
%         You must also select the variable algorithms from the following options:
%
%         algorithms={'active-set','interior-point','sqp','pso','ga',...
%                     'patternsearch'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear, close all
%parpool('local')
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOADING THE MEDITERRANEAN CORRIDOR
%%% MC structure containing the data for mediterranean corridor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('DATA/MC.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ( Model -> Scenarios 1, 2 and 3 ) NameModel={'LO', 'UP', '0'}
%    creditValue=50 euros and 23 g.tm by railway
%                          Environmental scenarios
%    149.7  g.tm by truck 'UP' (Scenario 1)
%     54    g.tm by truck 'LO' (Scenario 2)
%     0                    '0' (Scenario 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose a scenario

NameModel='0';

if strcmp(NameModel,'UP')
    N.eta=((149.7-23)/(10^6))*54.21;
elseif strcmp(NameModel,'LO')
    N.eta=((54-23)/(10^6))*54.21;
elseif strcmp(NameModel,'0')
    N.eta=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  LOAD EXISTING DATA TO ADD NEW RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['./RESULTS/experiment1_' NameModel '.mat'],'lambda_opt','f_opt',...
    'E_opt','N_opt','CPUtime','coste_igual','pr')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Units: kilometers, tonnes, and hours.
%% E structure associated with the events
%% N structure associated with network loading (solution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N.Tmax=24*365;                   % Planning period 24 hours * 365 days (1 year)
%  selection of a set of origin-destination pairs
N.nOD=1:length(MC.demand.L);     % All pairs origin-destination
%N.nOD= [140 260 604]            % a sample of o-d pairs
%N.nOD=1:50:600;                 % a sample of o-d pairs
N.Delta=1;                       % discretization = 1 train (size of package)
N.nTrainHour=6;                  % capacity parameter k_a Trains/hour per line
N.CommercialSpeed=53;            % km/h
N.Speed=100;                     % Speed in the running section (km/h)
N.FreightVolume=1230;            % Tm for a train of 750 meters
N.tau_a=G.Edges.Weight/N.Speed;  % travel time in the links


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          OPTIMIZACION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% auxiliary function needed to use a subset of o-d pairs
lambda=zeros(size(N.nOD ));
for i=1:length(N.nOD)
    indices(N.nOD(i))=i; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL 0:   PROPORTIONAL PRICES (results save in the index i=6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% definition of the objective function
i=6;                                         %  proportional approach
f=@(x) ObjectiveFunction(x,indices,N,MC,G);
expand_x=@(x) repmat(x,[1,length(N.nOD)])';
g=@(x) f(expand_x(x));                       % apply the single price to all o-d pairs

%% optimization
CPU_0=toc;
[lambda_opt{i},f_opt(i)] = fminbnd(g,0,0.5); % optimization of a single variable
CPUtime(i)=toc-CPU_0;                        % CPU time necessary to optimize
[Z,Z1,Z2,E_opt{i},N_opt{i}]=g(lambda_opt{i});% simulation at optimal prices

%% Plotting the revenue vs price
pr=0:0.01:0.5;                               % tabular data of prices
for j=1:length(pr)
    [coste_igual(j,1),coste_igual(j,2),coste_igual(j,3)]=g(pr(j));
end
plot(pr,coste_igual(:,1),pr,coste_igual(:,2),pr,coste_igual(:,3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL 1:   PATH-BASED PRICES (results saved in indices i=1:5, 7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = zeros(size(N.nOD ));                 % Constraints for prices. 
ub = 0.25*ones(size(N.nOD ));             % They are between 0 and 25% of fixed costs
lambda0=lambda_opt{6}*ones(size(N.nOD )); % initialization to proportional-based solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% selection of algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%algorithms={'active-set','interior-point','sqp','pso','ga','patternsearch'}
algorithms={'interior-point'};

for j=1:length(algorithms)

    if strcmp(algorithms{j},'active-set')
        i=1;
    elseif strcmp(algorithms{j},'interior-point')
        i=2;
    elseif strcmp(algorithms{j},'sqp')
        i=3;
    elseif strcmp(algorithms{j},'pso')
        i=4;
    elseif strcmp(algorithms{j},'ga')
        i=5;
    elseif strcmp(algorithms{j},'patternsearch')
        i=7;
    end

    CPU_0=toc;

    if any(strcmp(algorithms{j},{'active-set','interior-point','sqp'}))

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% EXACT METHODS 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        options = optimoptions('fmincon','Display','iter','Algorithm',algorithms{j},...
            'MaxFunctionEvaluations',4000);

        [lambda_opt{i},f_opt(i)] = fmincon(f,lambda0,[],[],[],[],lb,ub,[],options);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PSO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif  strcmp(algorithms{j},{'pso'})

        options = optimoptions('particleswarm','SwarmSize',60,'MaxIterations',500,...
            'Display', 'iter','UseParallel',true);

        [lambda_opt{i},f_opt(i)] = particleswarm(f,length(N.nOD) ,lb,ub,options);

        options1 = optimoptions('fmincon','Display','iter','Algorithm','active-set',...
            'MaxFunctionEvaluations',4000);

        [lambda_opt{i},f_opt(i)] = fmincon(f,lambda_opt{i},[],[],[],[],lb,ub,[],options1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GA 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp(algorithms{j},{'ga'})

        options = optimoptions('ga','UseParallel',true,'PopulationSize',60,'MaxGenerations',500,...
            'Display', 'iter');

        [lambda_opt{i},f_opt(i)] = ga(f,length(N.nOD),[],[],[],[],lb,ub,[],[],options);

        options1 = optimoptions('fmincon','Display','iter','Algorithm','active-set',...
            'MaxFunctionEvaluations',4000);

        [lambda_opt{i},f_opt(i)] = fmincon(f,lambda_opt{i},[],[],[],[],lb,ub,[],options1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% patternsearch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp(algorithms{j},{'patternsearch'})

        options = optimoptions('patternsearch','UseParallel',true,'Display','iter',...
            'MaxFunctionEvaluations',75000);

        [lambda_opt{i},f_opt(i)]=patternsearch(f,lambda0,[],[],[],[],lb,ub,[],options);

    end

    CPUtime(i)=toc-CPU_0;
     % simulation at optimal price
    [Z,Z1,Z2,E_opt{i},N_opt{i}]=ObjectiveFunction(lambda_opt{i},indices,N,MC,G);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['./RESULTS/experiment1_' NameModel '.mat'],'lambda_opt','f_opt','E_opt','N_opt','CPUtime','coste_igual','pr')

