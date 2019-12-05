% Script para gerar os resultados das simulações 
% das numerologias 1 e 3 utilizando as tecnicas de alocação
% Estatica - Dinamica para Máxima Vazão - e, AWM/MOM

%% Loading Results

% Static Allocation Num 1
SimData=load('staticAllocation_num1.mat');
EA1 = SimData.Sim.DataSNR;
EA2 = SimData.Sim.DataBPRB;  

% Static Allocation Num 3
SimData=load('staticAllocation_num3.mat');
EA3 = SimData.Sim.DataSNR;
EA4 = SimData.Sim.DataBPRB;

% Dynamic Allocation Num 1
DynamicData=load('dynamicMaxVazao_num1.mat');
DA1 = DynamicData.Dynamic.DataSNR;
DA2 = DynamicData.Dynamic.DataBPRB;

% Dynamic Allocation Num 3
DynamicData=load('dynamicMaxVazao_num3.mat');
DA3 = DynamicData.Dynamic.DataSNR;
DA4 = DynamicData.Dynamic.DataBPRB;

% AWM/MOM - Num 1
DynamicMOM=load('dynamic_mom_tese_num1.mat');
AWM1 = DynamicMOM.DynamicMOM.DataSNR;
AWM2 = DynamicMOM.DynamicMOM.DataBPRB;

% AWM/MOM - Num 3
DynamicMOM=load('dynamic_mom_tese_num3.mat');
AWM3 = DynamicMOM.DynamicMOM.DataSNR;
AWM4 = DynamicMOM.DynamicMOM.DataBPRB;


%% Gera graficos Alta Vazao
figure;
plot(EA1,EA2,'-.rs',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',6)
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;
hold on;

plot(DA1, DA2, '-.b*');

plot(AWM1, AWM2, '-.ko');

lgd0 = legend('Alocação Estática','Alocação Dinâmica/Máx. Vazão', 'Aloc. Dinâmica/AWM-MOM');
lgd0.FontSize = 10;
%% Gera graficos IoT
figure;
plot(EA3,EA4,'-.rs',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',6)
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;
hold on;
plot(DA3, DA4, '-.b*');

plot(AWM3, AWM4, '-.ko');

lgd = legend('Alocação Estática','Aloc. Dinâmica/Máx. Vazão', 'Aloc. Dinâmica/AWM-MOM');
lgd.FontSize = 10;