%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        MainExperiment1.m processes the experimental data 
%%%        corresponding to the section "Results I: Establishing 
%%%        Charges from the IM%s Perspective" in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         READING EXPERIMENTAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc,clear, close all
load('DATA/MC.mat'); % loading the case of study (Mediterranean Corridor)

% Three scenarios are considered, coded in the cell array Esc
% UP -> Scenario 1 in the paper
% LO -> Scenario 2 in the paper
% 0  -> Scenario 3 in the paper

Esc={'UP','LO','0'};  % Define the scenarios
Titulo={'Scenario 1','Scenario 2','Scenario 3'};
for escenario=1:length(Esc)

    NameModel=Esc{escenario};
    load(['./RESULTS/experiment1_' NameModel '.mat'],'lambda_opt','f_opt','E_opt','N_opt','CPUtime','coste_igual','pr')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      Table 5.  Computational costs   Path-based apporach 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tabla=[ [f_opt(1:5),f_opt(7),f_opt(6)]; [CPUtime(1:5), CPUtime(7), CPUtime(6)]];
    tabla=[round(tabla(1,:),2);
        round(tabla(2,:),1)];
    tabla(1,:)=-tabla(1,:);
    T=table(tabla(:,1),tabla(:,2),tabla(:,3),tabla(:,4),tabla(:,5),tabla(:,6),tabla(:,7),...
        'VariableNames',{'Active-set','Interior-point','SQP','PSO','GA','PatternSearch','Proportional'},...
        'RowNames',{'Z','CPU(s)'})
    table2latex(T, ['./RESULTS/tablaExperiment1_' NameModel '.tex'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      Figure 5: Objective function for the proportional approach
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    coste_igual(:,1)=-coste_igual(:,1);%
    ax=gca;
    ax.LineWidth=2;
    ax.FontSize=25;
    ax.Box='on';
    formato={'r*-','k^-','bo-'};
    for j=1:3
        hold on
        p=plot(pr',coste_igual(:,j),formato{j});
        p.LineWidth=2;
        xlabel('p');
        ytxt =char(8364);
        ylabel(['M ' num2str(ytxt)]);
        grid on
        title(Titulo{escenario})
    end
    legend({'Revenues (z)','Access fee','Environmental'})
    NameDibujo=['ProportionalCost' NameModel]
    ExportaGrafica(NameDibujo,ax)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         Table 6 Numerical results
    % %%        Comparision of PATHS-BASE COSTS vs PROPORTIONAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----
    i=6 % PROPORTIONAL
    E=E_opt{i};
    N=N_opt{i};
    tabla=procesar(E,N,MC,G);


    %---- PATH BASED
    [z,i_opt]=min([f_opt(1:5),Inf,f_opt(7)]) % BEST OF THE PATHS-BASED COSTS
    E=E_opt{i_opt};
    N=N_opt{i_opt};
    tabla=[procesar(E,N,MC,G);tabla];

    T=[round(tabla(1,:),2);
        round(tabla(2,:),2)];

    T=table(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),T(:,6),T(:,7),...
        'VariableNames',{'Access','CO2e','Rail cost','Delay cost','Speed Ave.','Tons (M)','Rail Share Market'},...
        'RowNames',{'Path-based','Proportional'})
    table2latex(T, ['./RESULTS/tablaCostExperiment1_' NameModel '.tex'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 6: Speed profiles vs scenarios (E,N PATH-BASED solution)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig=figure
    index_state_2=find(E.state==2);
    nbinds=200;
    subplot(1,2,1), hist(E.v(index_state_2),nbinds)
    axis([0,100,0,400])
    ax=gca;
    ax.LineWidth=2;
    ax.FontSize=18;
    ax.Box='on';
    ylabel('Number of trains')
    xlabel('Speed (km/h)')
    %hold on
    pd = fitdist(E.v(index_state_2)','Kernel','Kernel','epanechnikov')
    x_values = 0:1:100;
    y = pdf(pd,x_values);
    subplot(1,2,2), plot(x_values,y,'LineWidth',3)
    ax=gca;
    ax.LineWidth=2;
    ax.FontSize=18;
    xlabel('Speed (km/h)')
    ylabel('f.d.p')
    ax.Box='on';
    title(Titulo{escenario})
    NameDibujo=['SpeedProfile' NameModel]
    ExportaGrafica(NameDibujo,fig)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 7: Optimal flows for the scenarios (E,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig=figure
    geolimits([30 60],[-13 43])
    geobasemap none
    coeficiente=16616;
    WidthPropor=5*(N.nT_a/coeficiente)+2.5;

    color=colormap(copper);
    color=color(size(color,1):-1:1,:);
    min_color=min(WidthPropor);
    max_color=max(WidthPropor);

    for i=1:MC.link.n_links
        hold on
        a=G.findedge(MC.link.link_o{i},MC.link.link_d{i});
        o=find(strcmp(MC.node.node,MC.link.link_o{i})==1);
        d=find(strcmp(MC.node.node,MC.link.link_d{i})==1);
        P1=[MC.node.lat(o) MC.node.lat(d)];
        P2=[MC.node.ing(o) MC.node.ing(d)];
        %geoplot(P1,P2,'LineWidth',WidthPropor(a),'Color',[.6 0 0])
        id_color=floor(1+ ( (WidthPropor(a) - 2.5) /(6.0026-2.5)) *255);
        geoplot(P1,P2,'LineWidth',4,'Color',color(id_color,:))

    end
    hold on
    geoscatter(MC.centroid.lat,MC.centroid.ing,30,'MarkerFaceColor','k')
    %for i=1:length(MC.centroid.node)
        %text(MC.centroid.lat(i),MC.centroid.ing(i),MC.centroid.node{i},'VerticalAlignment','cap');
    %end
    geolimits([34 50],[-7 24])
    fontsize(fig, scale=1.2)
    title(Titulo{escenario})
    NameDibujo=['Flows' NameModel]
    ExportaGrafica(NameDibujo,fig)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save the plot in a file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig=DibujaSolucion(lambda_opt,Titulo{escenario})
    NameDibujo=['Solucion' NameModel]
    ExportaGrafica(NameDibujo,fig)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           A u x i l i a r   f u n c t i o n s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing of the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=procesar(E,N,MC,G)
index_state_1=find(E.state==1);
index_state_2=find(E.state==2);
index_state_0=find(E.state==0);
fprintf( 'Packages travelled = %i \n' , length(index_state_2))
fprintf( 'Packages in the newtork = %i \n' , length(index_state_1))
mTons=(N.FreightVolume*length(index_state_2))/10^6;
x(6)=mTons;
fprintf( 'Millons of tons by railway= %5.1f \n' , mTons)
ShareMarket=(mTons/MC.demand.totalD)*100;
x(7)=ShareMarket;
fprintf( 'Railway Share Market= %5.1f \n' , ShareMarket)
x(5)=(mean(E.v(index_state_2)));
fprintf(' Average speed = %5.1f \n ', x(5))

fprintf(' ---------- C O S T E S  ---------- \n ')
x(1)=(N.AccessCharges/10^6);
fprintf(' --- A C C E S O = %5.1f \n ',x(1))
x(4)=(N.cT/10^6);
fprintf(' ---  D E M O R A = %5.1f \n ',x(4))
x(3)=(N.Cr/10^6);
fprintf(' ---  T R A N S P O R T E = %5.1f \n ',x(3))
x(2)=(N.CO2e/10^6);
fprintf(' ---  E M I S S I O N S = %5.1f \n ',x(2))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export the generated plot to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportaGrafica(name,ax)
exportgraphics(ax,['./FIGURES/' name '.pdf'],'ContentType','vector')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig=DibujaSolucion(lambda_opt,Titulo)
proportional=lambda_opt{6}*ones([1,length(lambda_opt{1})]);
fig=figure
plot(lambda_opt{7},'r','LineWidth',2)
hold on
plot(proportional,'k.','LineWidth',2)
axis([0,630,0,0.25])
ax=gca;
ax.LineWidth=2;
ax.FontSize=18;
ax.Box='on';
xlabel('O-D pair')
ylabel('p')
title(Titulo)
lgd = legend({'Path-based .','Proportional'},'Location','southwest');
title(lgd,'Approach')

end
