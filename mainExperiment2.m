%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        MainExperiment2.m processes the experimental data 
%%%       corresponding to the section "Results II: Establishing the 
%%%       Pricing Method from the State%s Perspective" in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc,close all,clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Costs of externalities according to different studies in the literature
%%% Processing these values is obtained Table 7: Aggregated cost figures of
%%% the negative externalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

EXTER_road=[0.89      0.89      0.26     0.26     0.28    0.28     2.26     2.26     0.43   0.43;
               0      0.38      0.04     0.07     0       0.02     0,       9.21     0.01   1.27;
            0.05      8.36      0.05     8.36     0       2.37     0.24     0.24     0.05   0.89;
            0.31      0.31      0.31     0.31     0          0     0.13     0.13     0.32   0.32;
               0      2.39      0.07     0.07     0.02    0.02     0.01     0.29     0.07   0.26];

EXTER_train=[0.46     0.46      0.46     0.46     0.09    0.09     0        0        0.14   0.14;         
            0         0.13      0.02     0.02     0       0.01     0        1.3      0.03   0.03;
            0         0.07      0        0.09     0       0        0.01     0.01     0.05   0.05;
            0.05      0.05      0        0        0       0        0        0        0.03   0.03;
            0.02      1.13      0        0        0.01    0.01     0        0       0.01    0.03;];


T_road=[sum(EXTER_road(:,[1 3 5 7 9])') ; sum(EXTER_road(:,[2 4 6 8 10])')]';
E_road=mean([sum(EXTER_road(:,[1 3 5 7 9])') ; sum(EXTER_road(:,[2 4 6 8 10])')]');
T_train=[sum(EXTER_train(:,[1 3 5 7 9])') ; sum(EXTER_train(:,[2 4 6 8 10])')]';
E_train=mean([sum(EXTER_train(:,[1 3 5 7 9])') ; sum(EXTER_train(:,[2 4 6 8 10])')]');

T_ext=[T_road T_train];
aux=round([E_road,E_train],2);
T_ext=[T_ext;aux];
T_ext=table(T_ext(:,1),T_ext(:,2),T_ext(:,3),T_ext(:,4),...
    'VariableNames',{'Lower','Up','Lower_t','Up_t'},...
    'RowNames',{'1','2','3','4','5','mean'})

table2latex(T_ext, ['./RESULTS/tablaExperiment2_EXTERNA.tex']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing the Figure 9: Investment Plan in NPV and the annual invesment mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inversion=[1362, 8523, 30447 72415, (32498+1123)];
periodo=[3, 4, 5, 5, 15];
t=2014:(2014+sum(periodo)-1);
inv_anual=inversion./periodo;
tasa_descuento=2.5/100; % depreciation of the value of money
s0=find(t==2023); % reference year for the value of money
tasa_impuesto=0.32% Proportion of invested money recovered by the state through taxation
ia=[];
s=0;
for i=1:length(periodo)
    for j=1:periodo(i)
    ia=[ia inv_anual(i)*(1+tasa_descuento)^(s0-s) ];
    s=s+1;
    end
end
ia=ia*(1-tasa_impuesto);
InvAnualMedia=1/(length(ia)-1)*trapz(ia)

figure
ax=gca;
ax.LineWidth=2;
ax.Box='on';
p=plot(t,ia,'r*-')
p.LineWidth=2;
ax.FontSize=18;
    xlabel('Year');
    ytxt =char(8364);
    ylabel(['M ' num2str(ytxt)]);
    grid on
    title('Invesment Plan in the RFC6')
    exportgraphics(ax,['./FIGURES/Invesment.pdf'],'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Table 8: Benefits (in M e) vs scenarios
%%        Table 10: Benefit Cost Ratio
%%%       PATHS-BASE COSTS vs PROPORTIONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         READING EXPERIMENTAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('DATA/MC.mat'); % loading the case of study (Mediterranean Corridor)

% Three scenarios are considered, coded in the cell array Esc
% UP -> Scenario 1 in the paper
% LO -> Scenario 2 in the paper
% 0  -> Scenario 3 in the paper

Esc={'0','LO','UP'};% Define the scenarios
tabla=[];
Titulo={'Scenario 1','Scenario 2','Scenario 3'};

for escenario=1:length(Esc)

 NameModel=Esc{escenario};
 load(['./RESULTS/experiment1_' NameModel '.mat'],'lambda_opt','f_opt','E_opt','N_opt','CPUtime','coste_igual','pr')

%----
i=6; % PROPORTIONAL
E=E_opt{i};
N=N_opt{i};
%tabla=procesar(E,N,MC,G);
tabla=[procesar(E,N,MC,G);tabla];

%---- PATH BASED
[z,i_opt]=min([f_opt(1:5),Inf,f_opt(7)]); % BEST OF THE PATHS-BASED COSTS
E=E_opt{i_opt};
N=N_opt{i_opt};
tabla=[procesar(E,N,MC,G);tabla];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing the benefits 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(1+tasa_descuento).^13.
Externalidad=E_road-  E_train;
%a valor social
a=(-1+(0.385/0.045)).*tabla(:,3).*(1-tasa_impuesto).*((tabla(:,7)-11.41)/100);
b1= (1+tasa_descuento).^13.*((aux(1)- aux(3))/100).*...
(tabla(:,2)-500*MC.demand.totalD.*0.1141);
b2= (1+tasa_descuento).^13.*((aux(2)- aux(4))/100).*...
(tabla(:,2)-500*MC.demand.totalD.*0.1141);

Beneficio1=tabla(:,1)+b1-tabla(:,4)+a;
Beneficio2=tabla(:,1)+b2-tabla(:,4)+a;

T=[round(tabla(:,1),2) round(b1,2) round(b2,2) round(-tabla(:,4),2)  round(a,2)...
  round(Beneficio1,2) round(Beneficio2,2)  round(100*(Beneficio1/InvAnualMedia),2)... 
  round(100*(Beneficio2/InvAnualMedia),2)];

T1=table(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),T(:,6),T(:,7),...
    'VariableNames',{'Access','Externalities_LO','Externalities_UP', ...
    'FOC''s benefit','Social'' benefit','Total Benefit_LO','Total Benefit_UP'},...
   'RowNames',{'1Path-based','1Proportional','2Path-based','2Proportional','3Path-based','3Proportional'});

T2=table(T(:,8),T(:,9),...
    'VariableNames',{'BCR_LO','BCR_UP'},...
   'RowNames',{'1Path-based','1Proportional','2Path-based','2Proportional','3Path-based','3Proportional'});

table2latex(T1, ['./RESULTS/tablaExperiment2Benefit.tex']);
table2latex(T2, ['./RESULTS/tablaExperiment2BCR.tex']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           A u x i l i a r   f u n c t i o n s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=procesar(E,N,MC,G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index_state_1=find(E.state==1);
index_state_2=find(E.state==2);
index_state_0=find(E.state==0);

mTons=(N.FreightVolume*length(index_state_2))/10^6;
x(6)=mTons;
%fprintf( 'Millons of tons by railway= %5.1f \n' , mTons)
ShareMarket=(mTons/MC.demand.totalD)*100;
x(7)=ShareMarket;
%fprintf( 'Railway Share Market= %5.1f \n' , ShareMarket)
x(5)=(mean(E.v(index_state_2)));
%fprintf(' Average speed = %5.1f \n ', x(5))

%fprintf(' ---------- C O S T E S  ---------- \n ')
x(1)=(N.AccessCharges/10^6);
%fprintf(' --- A C C E S O = %5.1f \n ',x(1))
x(4)=(N.cT/10^6);
%fprintf(' ---  D E M O R A = %5.1f \n ',x(4))
x(3)=(N.Cr/10^6);
%fprintf(' ---  T R A N S P O R T E = %5.1f \n ',x(3))
x(2)=((N.tnkm)/10^6);
%fprintf(' ---  E M I S S I O N S = %5.1f \n ',x(2))

end


