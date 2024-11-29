%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         O B J E C T I V E      F U N C T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,Z1,Z2,E,N]=ObjectiveFunction(lambda,indices,N,MC,G)
N.Cr=0;                               % fixed cost (euros)
N.CO2e=0;                             % value of the externality cost
N.tnkm=0;                             % ton X km in the rail network
N.AccessCharges=0;                    % access cost (euros)
N.cT=0;                               % delay cost by (euros)
N.nT_a=zeros(size(N.tau_a));          % number of package in the arc
N.n_package_i=zeros([1,max(N.nOD)]);  % number of packages in each OD
AccessCostPricing=@(T,r) lambda(indices(r));
[E,N]=Simulation(AccessCostPricing,N,MC,G);
Z1=N.AccessCharges / (10^(6));        % Revenue from railway network access charges (First term of objetive function)
Z2=N.CO2e / 10^(6);                   % Externality term (Second term of the objective function)
Z=-(Z1+Z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      S I M U L A T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E,N]=Simulation(AccessCostPricing,N,MC,G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTIALIZATION OF EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E.np=0;                                    % counter  of packages;
E.ne=0;                                    % counter of events
E.list=[];                                 % list of events
T=0;                                       % simulation clock
for a=1:G.numedges
    for s=1:2
        E.q{a,s}=[];                       % Vertical queues of the links
    end
end

for r=N.nOD
    E.ne=E.ne+1;                           % update event counter
    j=E.ne;                                % new event
    E.list(end+1)=j;                       % add event to list of events
    E.T(j)=0;                              % instant of the event
    E.a(j)=1;                              %  arc (link) of the instant
    E.s(j)=1;                              % vertex of the event
    E.t0(j)=0;                             %  entrance time
    E.r(j)=r;                              % only a path,, i.e  r=omega
    E.state(j)=0;                          % state of the event =0 (to entrance to the network) 
                                           %                    =1 (in the network)
    N.tau(r)=MC.demand.L(r)/N.CommercialSpeed; %current  travel time for path r
    N.bar_tau(r)=MC.demand.L(r)/N.CommercialSpeed; %Reference  travel time for path r
    N.kappa(r)=1/N.FreightVolume;      %  rate to convert freight flow to train flow;
    N.v(r)=N.CommercialSpeed;            % current speed for path r ( reference speed for path r)
    N.L(r)=MC.demand.L(r);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATION MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ( T<=N.Tmax &  not(length(E.list)==0) )

    [T,i_indx]=min(E.T(E.list));            % update the simulation clock
    i=E.list(i_indx);

    if E.state(i)==0                                    % STATE0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new ficiticious event j
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        E.ne=E.ne+1;
        j=E.ne;
        E.list(end+1)=j ;
        E.r(j)=E.r(i);
        E.state(j)=0;
        E.a(j)= 1;                                      %  initial link
        E.s(j)= 1;                                      %  initial vertex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Loading package i onto the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        E.np=E.np+1;
        E.state(i)=1;
        o=NumNode(MC.node.iso2{find(strcmp(MC.node.node,MC.demand.o{E.r(i)})==1)});
        d=NumNode(MC.node.iso2{find(strcmp(MC.node.node,MC.demand.d{E.r(i)})==1)});
        pr=AccessCostPricing(E.T(i),E.r(i));

        % Railway cost
        [A_r,lambda_r]=CostRailway(E.r(i),MC.demand.LogitModel, N,pr);
        
        N.AccessCharges=N.AccessCharges+lambda_r*N.tau(E.r(i))*N.L(E.r(i))*N.FreightVolume;
        N.cT=N.cT+A_r*N.tau(E.r(i))*N.L(E.r(i))*N.FreightVolume;
        N.Cr=N.Cr+N.L(E.r(i))*MC.demand.LogitModel.CostRailKilometer*N.FreightVolume;
        N.CO2e=N.CO2e+N.eta*N.L(E.r(i))*N.FreightVolume;
        N.tnkm=N.tnkm+N.L(E.r(i))*N.FreightVolume;

        % Share Market
        [RailwayShareMarket] = Logit_r(o,d,MC.demand.LogitModel,A_r,lambda_r,N.tau(E.r(i)));

        IncT=1/(RailwayShareMarket*MC.demand.D(E.r(i))*N.kappa(E.r(i)));

        E.T(i)=E.T(i)+IncT;
        E.t0(i)=E.T(i);
        E.T(j)=E.T(i);
    end

    if  E.a(i) <MC.demand.n_links(E.r(i))

        % queue section

        a=MC.demand.link_path{E.r(i)}(E.a(i));              % id link
        s=MC.demand.node_path{E.r(i)}(E.s(i));              % name_node
        id_s=strcmp(MC.demand.Edges.EndNodes(a,1),s)+1;     % direction of the link 1 / 2

        Time_queue=CapacityDynamic(MC.link.link_capacity(a),N.Delta,E.T(i),N.nTrainHour);
      
        for j=E.q{a,id_s}                                  % for all package in the queue
            if E.T(j)<(E.T(i)+Time_queue)
                E.T(j)=E.T(i)+Time_queue;
            end
        end

        id_i=find(E.q{a,id_s}==i);
        E.q{a,id_s}(id_i)=[];                             % Delete i from the queue
        % running section

        E.T(i)=E.T(i)+N.tau_a(a); % Travel time in the running section
        N.nT_a(a)=N.nT_a(a)+1;    % pass a package trough the arc a

        % advance the package i in the network
        E.a(i)=E.a(i)+1;                                  %  next link
        E.s(i)=E.s(i)+1;                                  %  next vertex
        a_mas=MC.demand.link_path{E.r(i)}(E.a(i));        % id link
        s_mas=MC.demand.node_path{E.r(i)}(E.s(i));        % name_node
        id_s_mas=strcmp(MC.demand.Edges.EndNodes(a_mas,1),s_mas)+1; % direction of the link (1/2)
        E.q{a_mas,id_s_mas}=[E.q{a_mas,id_s_mas} i]; %    queue i

    else
        E.state(i)=2;                 % arrive the package to destination
        id_i=find(E.list==i);
        E.list(id_i)=[];
        E.tau(i)=E.T(i)- E.t0(i);
        E.v(i)=MC.demand.L(E.r(i))/E.tau(i);
        % average speed
        % update the travel time in path r
        N.n_package_i(E.r(i))=N.n_package_i(E.r(i))+1;
        p1=1/N.n_package_i(E.r(i));
        N.v(E.r(i))= (1-p1)*N.v(E.r(i))+p1*E.v(i);
        N.tau(E.r(i))=(1-p1)* N.tau(E.r(i))+p1* E.tau(i);
    end
end
end