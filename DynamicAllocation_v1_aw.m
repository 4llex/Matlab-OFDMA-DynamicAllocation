%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

%% Define Numerology
Numerology = 3;

if (Numerology == 1)
     N = 6336;
     sc_per_rb = 48;
     RE = 720;
else
     N = 1584;
     sc_per_rb = 12;
     RE = 720;
end

%%
TargetSer = 1e-3;                           %% SER Alvo
SNR = 6:2:30;                               %% XXX
%N = 6336;                                  %% Numero de Subportadoras
b = zeros(1,N);                             %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));          %% Total de bits em um simbolo
bits_per_rb = zeros(1,length(SNR));         %% qtd media de Bits por RB 
quantizar = 'yes';                          %% 
RB = 132;                                   %% qtd de RB
%sc_per_rb = 48;                            %% SubCarriers per RB, depends numerology    
nusers = 3;
%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao
% LTE EVA CHANNEL
freq_sample = 23.76e6;     %N*15e3; %30.72e6; sample rate do LTE
EVA_SR3072_Delay           =[0 30 150 310 370 710 1090 1730 2510].*1e-9;
EVA_SR3072_PowerdB_Gain    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; %  20*log10(0.39)= -8.1787 => -8.1787 dB -> Voltage-ratio = 0.398107

chan_EVA = rayleighchan((1/(freq_sample)),0,EVA_SR3072_Delay,EVA_SR3072_PowerdB_Gain);        
impulse= [1; zeros(N - 1,1)];  


H    = ones(nusers,RB);
mask = zeros(nusers,RB);
capacity = zeros(nusers,RB);


num_itr = 200;
for i=1:length(SNR)
    i
    j=0;
    while j<num_itr 
        
        % Gera o canal randomico para cada user
        for user=1:nusers
            h = filter(chan_EVA, impulse)';
            Hf = fft(h,N);
            % Calcula Resposta em frequencia média para os 132 RB's
            H(user,:) = rb_h_media(Hf, sc_per_rb);
        end
        
        SNRLIN = 10^(SNR(i)/10);
        P  = 20;
        Pu = P/nusers;
        
        for user=1:nusers
            mask(user,:) = ( abs(H(user,:))== max(abs(H)) );
            [~,~, capacity(user,:) ] = fcn_waterfilling(Pu, P/(SNRLIN*RB), Gamma, H(user,:), mask(user,:) );
        end
%       sum(mask(:))

   
        b = sum(capacity);

        % Quantização
        if strcmp(quantizar,'yes')
            b(b<2) = 0;
            b((b>2)&(b<4)) = 2;
            b((b>4)&(b<6)) = 4;
            b((b>6)&(b<8)) = 6;
            b(b>8) = 8;
        end

%         figure(2);
%         cla;
%         plot(subPower1); hold on; plot(subPower2); plot(subPower3);
%         drawnow();

%         figure(3);
%         cla;
%         plot(abs(H1)); hold on; plot(abs(H2)); plot(abs(H3));
%         drawnow();
%         pause(1);
        
        Total_bits(i) = Total_bits(i) + sum(b);
        %bm(j) = b;
        
        j = j+1;
    end  
    
    Total_bits(i) = Total_bits(i)/num_itr;
    
    bits_per_rb(i) = (Total_bits(i)/RB)*RE; 
end

%% Loading File
if (Numerology == 1)
     SimData=load('staticAllocation_num1.mat');
     D1 = SimData.Sim.DataSNR;
     D2 = SimData.Sim.DataBPRB;  
else     
     SimData=load('staticAllocation_num3.mat');
     D1 = SimData.Sim.DataSNR;
     D2 = SimData.Sim.DataBPRB;
end  

%% Gera graficos de Bits/SNR
figure;
plot(SNR, bits_per_rb, '-o');
%title('Alocação de Recursos em sistema de multiplo acesso Ortogonal');
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;

hold on;
plot(D1, D2, '-or');
legend('Alocação Dinâmica','Alocação Estática')