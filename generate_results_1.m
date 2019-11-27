% Script para gerar os resultados das simulações 
% das numerologias 1 e 3 utilizando as tecnicas de alocação
% Estatica Vs Dinamica para Máxima Vazão (para numerologias 1 e 3)

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


%% Gera graficos de Bits/SNR
figure;
plot(EA1, EA2,'-.rd');
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;

hold on;
plot(EA3,EA4,'-.rs',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',6)
%plot(EA3, EA4, ':rs');

hold on;
plot(DA1, DA2, '-.bo');

hold on;
plot(DA3,DA4,'-.bv',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',6)
%plot(DA3, DA4, ':bv');

lgd = legend('Alocação Estática - Alta Vazão','Alocação Estática - IoT', 'Alocação Dinâmica - Alta Vazão', 'Alocação Dinâmica - IoT');
lgd.FontSize = 10;