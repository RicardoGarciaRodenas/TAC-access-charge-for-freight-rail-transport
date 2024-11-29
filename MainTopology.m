%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the case study for the Mediterranean Corridor RFC6.
% After execution, it generates a MATLAB file MC.mat that stores the structure
% (data type) MC and the graph G containing all the necessary information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure found at: .node
% Field: .node.node - Type: cell
% Field: .node.iso2 - Type: cell
% Field: .node.mainNode - Type: double
% Field: .node.Port - Type: double
% Field: .node.Terminals - Type: double
% Field: .node.RFC - Type: double
% Field: .node.lat - Type: double
% Field: .node.ing - Type: double

% Structure found at: .link
% Field: .link.link_o - Type: cell
% Field: .link.link_d - Type: cell
% Field: .link.link_capacity - Type: double
% Field: .link.link_dist - Type: double
% Field: .link.n_links - Type: double


% Structure found at: .centroid
% Field: .centroid.node - Type: cell
% Field: .centroid.iso2 - Type: cell
% Field: .centroid.Port - Type: cell
% Field: .centroid.TotalTrafficPort - Type: cell
% Field: .centroid.Terminals - Type: cell
% Field: .centroid.RFC - Type: cell
% Field: .centroid.lat - Type: double
% Field: .centroid.ing - Type: double
% Structure found at: .demand

% Field: .demand.flow - Type: double
% Structure found at: .demand.LogitModel
% Field: .demand.LogitModel.betaRail - Type: double
% Field: .demand.LogitModel.betaTruck - Type: double
% Field: .demand.LogitModel.alfa - Type: double
% Field: .demand.LogitModel.CostTruck - Type: double
% Field: .demand.LogitModel.CostRailKilometer - Type: double
% Field: .demand.LogitModel.CostRailTime - Type: double
% Field: .demand.totalD - Type: double
% Field: .demand.L - Type: double
% Field: .demand.o - Type: cell
% Field: .demand.d - Type: cell
% Field: .demand.node_path - Type: cell
% Field: .demand.link_path - Type: cell
% Field: .demand.n_links - Type: double
% Field: .demand.D - Type: double
% Structure found at: .demand.Edges
% Field: .demand.Edges.EndNodes - Type: cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MEDITERRANEAN RAIL FREIGHT CORRIDOR
%% Spain, France, Italy, Eslovakia, Croacia, Hungry,
% https://www.rfi.it/en/In-Europe/Freight-Corridors/Mediterranean-Freight-Corridor.html
% Routing: Almería-Valencia/Algeciras/Madrid-Zaragoza/Barcelona-Marseille-Lyon-Turin-Milan-Verona -Padua/Venice-Trieste/Capodistria-Ljubljana-Budapest Ljubljana/Fiume -Zagabria-Budapest-Zahony (Hungarian-Ukranian border);
% Members: ADIF (Spain), Línea Figueras Perpignan (Spain-France), SNCF Réseau (France), OcVia (France), RFI (Italy), SŽ - Infrastruktura (Slovenia), HŽ Infrastruktura (Croatia), MÁV (Hungary); VPE (Hungary), e .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear,clc,close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOADING RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DATA/worldcities.mat','data') % cities of countries of the mediterranea corridor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESSING THE NODES (yards and railway junctions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sheet=1; % nodes
[numNodes,txtNodes,rawNodes] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
[n_nodes,n_atributes]=size(numNodes);
for j=1:n_nodes
    MC.node.node{j}=txtNodes{j+1,1};
    MC.node.iso2{j}=txtNodes{j+1,8};    % country
    MC.node.mainNode(j)=numNodes(j,1);  % primary / secondary nodes
    MC.node.Port(j)=numNodes(j,2);
    MC.node.Terminals(j)=numNodes(j,3);
    MC.node.RFC(j)=numNodes(j,4);
    if isnan(numNodes(j,5))
        indx=find(strcmp(data.city, MC.node.node{j}));
        MC.node.lat(j)=data.lat(indx);
        MC.node.ing(j)=data.ing(indx);
    else
        MC.node.lat(j)=numNodes(j,5);   % latitud
        MC.node.ing(j)=numNodes(j,6);   % longitud
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESSING THE LINKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sheet=2; % links
[numLinks,txtLinks,rawLinks] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
[n_links,n_atributes]=size(txtLinks);
% processing of the vertices of the links
for j=1:n_links
    MC.link.link_o{j}=txtLinks{j,1};
    MC.link.link_d{j}=txtLinks{j,2};
    MC.link.link_capacity(j)=rawLinks{j,3};
end

%%%%%%%%%%%%%%%%%%
%% computing the length of the links
%%%%%%%%%%%%%%%%%%

wgs84 = wgs84Ellipsoid;
D=zeros([1,6]); %length of the link (kilometer)
for i=1:size(MC.link.link_o,2)
    o=find(strcmp(MC.node.node,MC.link.link_o{i})==1);
    d=find(strcmp(MC.node.node,MC.link.link_d{i})==1);
    P1=[MC.node.lat(o) MC.node.lat(d)];
    P2=[MC.node.ing(o) MC.node.ing(d)];
    MC.link.link_dist(i)=distance(MC.node.lat(o),MC.node.ing(o),MC.node.lat(d),MC.node.ing(d),wgs84)/1000;
    A=1;
    if not(isnan(numLinks(i,1)))
        A=2;
    end
    D(NumNode(MC.node.iso2{o}))=A*MC.link.link_dist(i) *numLinks(i,2)+D(NumNode(MC.node.iso2{o}));
    D(NumNode(MC.node.iso2{d}))=A*MC.link.link_dist(i) *numLinks(i,3)+D(NumNode(MC.node.iso2{d}));

end

D_ref=[3015, 1515 636 457 375 1143]; % total real length of the Mediterranean RFC (principal route)
corretion_rate= D_ref./D;

%%%%%%%%%%%
%% correction of the length
%%%%%%%%%%%

for i=1:size(MC.link.link_o,2)
    o=find(strcmp(MC.node.node,MC.link.link_o{i})==1);
    d=find(strcmp(MC.node.node,MC.link.link_d{i})==1);
    MC.link.link_dist(i)=MC.link.link_dist(i) *(numLinks(i,2)* corretion_rate(NumNode(MC.node.iso2{o}))+...
        numLinks(i,3)* corretion_rate(NumNode(MC.node.iso2{d})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTING THE CENTROIDS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sheet=3; % centroids
[numCentroids,txtCentroids,data_centroids] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
[n_centroids,n_atributes_centroid]=size(txtCentroids);

for j=1:(n_centroids-1)
    MC.centroid.node{j}= data_centroids{j+1,1}; %node
    MC.centroid.iso2{j}=data_centroids{j+1,8}; % country
    MC.centroid.Port(j)=data_centroids(j+1,4); % port?
    MC.centroid.TotalTrafficPort(j)=data_centroids(j+1,5); % traffic flow for port
    MC.centroid.Terminals(j)=data_centroids(j+1,6); % number of terminals
    MC.centroid.RFC{j}=data_centroids{j+1,7}; % zone geograffic to connect
    o=find(strcmp(MC.node.node,MC.centroid.node{j})==1);
    MC.centroid.lat(j)=MC.node.lat(o);
    MC.centroid.ing(j)=MC.node.ing(o);
    try
        indx=find(strcmp(data.city, txtCentroids{j+1}));
        Pob(j)=data.population(indx);

    catch
        Pob(j)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTING THE O-D FLOWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%%% projection matrix to dissagregate the flows
%%%%%%%%%%%

mTons=zeros([n_centroids-1,n_centroids-1]); % demand / year (millons of tons)
P_pob=ProjectionMatrix(Pob);
P_traff =ProjectionMatrix(cell2mat(MC.centroid.TotalTrafficPort));
P_term=ProjectionMatrix(cell2mat(MC.centroid.Terminals));
P_traff_tem=0.5*P_term+0.5*P_traff;

%%%%%%%%%%%
%%%  dissagregation of the flows. Countries to nodes
%%%%%%%%%%%

sheet=4; % OD target scenario 3
[numOD,txtOD,data_OD] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
[n_zones,n_atributes_zones]=size(txtOD);
sheet=5; % Railway Share Market
[numShareMarket,txtShareMarket,data_ShareMarket] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
% total freight flow for all modes of transport
TotalFlow=cell2mat(data_OD(2:end,2:end))./cell2mat(data_ShareMarket(2:end,2:end));
[i_max,j_max]=size(TotalFlow);
for i=1:i_max
    o=data_OD{i+1,1};
    for j=1:j_max
        d=data_OD{1,j+1};
        Flow_od=TotalFlow(i,j);
        if not(isnan(Flow_od))
            if o=='SE'
                indx_o=find(strcmp(MC.centroid.RFC, o));
            else
                indx_o=find(strcmp(MC.centroid.iso2, o));
            end
            if  any([strcmp(d,'SE' ),  strcmp(d,'NE' ) ,strcmp(d,'W' )])
                indx_d=find(strcmp(MC.centroid.RFC, d));
            else
                indx_d=find(strcmp(MC.centroid.iso2, d));
            end
            % disaggregation for zone by zone
            if  any([strcmp(o,'SE' ), strcmp(d,'SE' ),  strcmp(d,'NE' ) ,strcmp(d,'W' )])
                mTons=mTons+Dissagregation(indx_o, indx_d,P_traff_tem,Flow_od);
            else
                mTons=mTons+Dissagregation(indx_o, indx_d,P_pob,Flow_od);
            end
        end
    end
end

mTons=(0.5*mTons+0.5*mTons')/1000;
MC.demand.flow=mTons;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALIBRATION OF THE LOGIT MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sheet=5; % Railway Share Market Scenario 3 (best case)
[numShareMarket,txtShareMarket,data_ShareMarket3] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
sheet=7; % Railway Share Market Scenario 2 (worst case)
[numShareMarket,txtShareMarket,data_ShareMarket2] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
% parameters
Count_od=0;
CostTruck=0.385;    % cost per kilometer and ton by road
CostRail=0.045;     % cost per kilometer and ton by railway
IncCostEscenario=0.18; %increment of cost in the mode road (difference between scenario 3 and 2)

for i=1:(i_max-1)
    o=data_OD{i+1,1};
    for j=1:(j_max-3)
        d=data_OD{1,j+1};
        Flow_od=TotalFlow(i,j);
        if not(isnan(Flow_od))
            % scenario 3
            Count_od=Count_od+1;
            y(Count_od)=log(1/data_ShareMarket3{i+1,j+1}-1);
            X(Count_od,i)=1;
            X(Count_od,j+1)=1;
            betaRoad(Count_od)=CostTruck*(1+IncCostEscenario);
            betaRailway(Count_od)=-CostRail;
            interceptationRail(Count_od)=-1;
            % scenario 2
            Count_od=Count_od+1;
            y(Count_od)=log(1/data_ShareMarket2{i+1,j+1}-1);
            X(Count_od,i)=1;
            X(Count_od,j+1)=1;
            betaRoad(Count_od)=CostTruck;
            betaRailway(Count_od)=-CostRail;
        end
    end
end
tbl = table(betaRoad',X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),y', ...
    'VariableNames',{'betaRoad','ES','FR','IT',	'SL',	'HR',	'HU','OddRatio'});
mdl = fitlm(tbl,'OddRatio~betaRoad+ES+FR+IT+SL+HR+HU')
%plot(mdl)
parameters=mdl.Coefficients.Estimate;
parameters(1)=-mdl.Coefficients.Estimate(1)/CostRail;

MC.demand.LogitModel.betaRail=parameters(1);
MC.demand.LogitModel.betaTruck=parameters(2);
MC.demand.LogitModel.alfa=parameters(3:end);
MC.demand.LogitModel.CostTruck=CostTruck;
MC.demand.LogitModel.CostRailKilometer=CostRail;
MC.demand.LogitModel.CostRailTime=2.23 % euros per Tm and hour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDING THE GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sheet=8; % conncectors centroids-node
[numConnectors,txtConnectors,dataConnectors] = xlsread('./DATA/TopologyMediterraneanCorridor.xlsx',sheet);
[n_Connectors,n_atributes_Connectos]=size(txtConnectors);

NodeWhole=MC.node.node;
MC.link.n_links=size(MC.link.link_o,2);          % number of links without connectors
%Link_o_Whole=MC.link.link_o;
%Link_d_Whole=MC.link.link_d;
%weights=MC.link.link_dist;
for j=2:n_Connectors
    NodeWhole{end+1}=txtConnectors{j,2};
    MC.node.node{end+1}=txtConnectors{j,2};
    MC.link.link_o{end+1}=txtConnectors{j,1};
    MC.link.link_d{end+1}=txtConnectors{j,2};
    MC.link.link_capacity(end+1)=Inf;
    MC.link.link_dist(end+1)=0;
end
NodeTable=table(NodeWhole','VariableNames',{'Name'});
G = graph(MC.link.link_o, MC.link.link_d,MC.link.link_dist,NodeTable);
%%% indices physical arcs within links

%%% reordering the capacity of the arcs
capacity_aux=MC.link.link_capacity;
for i=1:length(MC.link.link_capacity)

    a=G.findedge(MC.link.link_o{i},MC.link.link_d{i});
    MC.link.link_capacity(a)=capacity_aux(i);
 
end


%MC.centroid.node
n_od=0;                                              
MC.demand.totalD=0;                         % total demand in the corredor
for i=1:length(MC.centroid.node)
    for j=1:length(MC.centroid.node)
        % Delete o-d pair with lesser than a train per year
        if ne(i,j) & (MC.demand.flow(i,j)>0) & (0.05*MC.demand.flow(i,j)*10^(6) >1230)
            n_od=n_od+1;
            o=['o_' MC.centroid.node{i}];
            d=['d_' MC.centroid.node{j}];
            [path,d,link_path] = shortestpath(G,o,d,'Method','positive');
            MC.demand.L(n_od)=d;
            MC.demand.o{n_od}=MC.centroid.node{i};
            MC.demand.d{n_od}=MC.centroid.node{j};
            MC.demand.node_path{n_od}=path;
            MC.demand.link_path{n_od}=link_path;
            MC.demand.n_links(n_od)=length(link_path);
            MC.demand.D(n_od)=(MC.demand.flow(i,j)*10^(6))/(365*24); % Tm per hour
            MC.demand.totalD= MC.demand.totalD+MC.demand.D(n_od);
        end
    end
end
MC.demand.totalD=(MC.demand.totalD*365*24)/(10^6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTTING THE MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%%%  PLOTTING THE CORRIDOR
%%%%%%%%%%%

figure(2)
geobasemap streets
% other options
%geobasemap topographic (nice)
%geobasemap none
%geobasemap grayland
hold on
for i=1:MC.link.n_links
    hold on
    o=find(strcmp(MC.node.node,MC.link.link_o{i})==1);
    d=find(strcmp(MC.node.node,MC.link.link_d{i})==1);
    P1=[MC.node.lat(o) MC.node.lat(d)];
    P2=[MC.node.ing(o) MC.node.ing(d)];
    geoplot(P1,P2,'LineWidth',2,'Color',[.6 0 0])
end
geolimits([35 55],[-10 25])

%%%%%%%%%%%
%%%  PLOTTING THE FLOW
%%%%%%%%%%%
figure(3)
geolimits([30 60],[-13 43])
geobasemap none
hold on
A=sum(mTons);
geoscatter(MC.centroid.lat,MC.centroid.ing,A*24,'MarkerFaceColor',[0.9290,0.6940,0.1250])
for i=1:length(MC.centroid.node)
text(MC.centroid.lat(i),MC.centroid.ing(i),MC.centroid.node{i},'VerticalAlignment','cap');
end

for i=1:MC.link.n_links
    hold on
    o=find(strcmp(MC.node.node,MC.link.link_o{i})==1);
    d=find(strcmp(MC.node.node,MC.link.link_d{i})==1);
    P1=[MC.node.lat(o) MC.node.lat(d)];
    P2=[MC.node.ing(o) MC.node.ing(d)];
    geoplot(P1,P2,'LineWidth',2,'Color',[.6 0 0])
end

geolimits([35 55],[-10 25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE  THE  MEDITERRANEAN CORRIDOR DATA (MC) AND THE GRAPH (G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MC.demand.Edges.EndNodes=G.Edges.EndNodes;

save('./DATA/MC.mat','MC','G')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               A U X I L I A R     F U N C T I O N S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=Dissagregation(o,d, P,F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D Demand
% o indices associated with origin
% d  indices associated with destination
% P projection matrix
D=zeros(size(P));
D(o,d)=P(o,d);
Total=sum(D,'all');
D=D./Total;
D=D.*F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P=ProjectionMatrix(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=x'*x;
P=P.*(ones(size(P))-eye(size(P)));
Total=sum(P,'all');
P=P./Total;
end
