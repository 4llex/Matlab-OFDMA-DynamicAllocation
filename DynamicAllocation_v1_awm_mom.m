%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

%% Water Filing Modificado para MOM, maxima ondem de modulação: 
%  cada subportadora sobressalente é alocada para o usuario que pode
%  atingir a maior quantidade de bits!

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
SNR = 10:2:40;                               %% XXX
%N = 6336;                                  %% Numero de Subportadoras
b = zeros(1,N);                             %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));          %% Total de bits em um simbolo
bits_per_rb = zeros(1,length(SNR));         %% qtd media de Bits por RB 
%quantizar = 'yes';                          %% 
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
%mask = zeros(nusers,RB);
capacity = zeros(nusers,RB);

% new variable for AWM
mask = ones(nusers,RB); % mask para WF em todas as portadoras, tudo em '1'
priority_user = zeros(1,nusers);
bmax = zeros(1,nusers);
%real_capacity = zeros(nusers,RB);
%test = [];

num_itr = 5000;
for i=1:length(SNR)
    i
    j=0;
    
    while j<num_itr 
        
        %bmin = [100, 100, 100];
        bmin = [50, 50, 50];
        
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
        
        % Distribuição de potencia utilizando WF, para cada user em todo o
        % espectro OFDM
        for user=1:nusers
            [~,~, capacity(user,:) ] = fcn_waterfilling(Pu, P/(SNRLIN*RB), Gamma, H(user,:), mask(user,:) ); % a mask é tudo '1'!
            bmax(user) = sum(capacity(user,:));
            capacity(user,:) = quantization(capacity(user,:));
        end
        
        % Gettting priority users
        for user=1:nusers
            [~,index] = max(bmax);
            priority_user(user) = index;
            bmax(index) = -1;
        end
        
        
        priority_user;
        alloc_vec = zeros(1, RB);
        alloc_user = zeros(1, RB);
        real_capacity = zeros(nusers,RB);
        
        while (sum(bmin<=0) ~= nusers)
            
                if(sum(alloc_vec)==132)
                    break;
                else
                    for ii=1:nusers
                        if (bmin(priority_user(ii))>0)
                           [value,index] = max(capacity(priority_user(ii),:));
                           real_capacity(priority_user(ii),index) = value; % value é ordem de modulação do um RB!
                           capacity(:,index) = -1;
                           alloc_vec(index) = 1;
                           alloc_user(index) = ii;
                           %test = [test,index];
                           bmin(priority_user(ii)) = bmin(priority_user(ii)) - value;
                        end
                    end % end do for
                end
        end % end while 
        
        % Verifica se há portadoras sobressalentes e aloca cada uma para o 
        % usuario que pode transmitir a maior taxa de bits!
        mask2 = zeros(nusers,RB); % mask para melhor user por portadora
        if (sum(alloc_vec)~=132)
            
            idx_sobressalentes = find(~alloc_vec); % retorna index das sc sobressalentes!
            %x= find(alloc_vec);
            for user=1:nusers
                mask2(user,idx_sobressalentes) = ( abs(H(user,idx_sobressalentes))== max(abs(H(:,idx_sobressalentes))) ); % mask é 1 onde o user pode transmitir melhor
            end
            %y = find(~(sum(mask2)));
            
            for user=1:nusers % Obtem idx da maskara de cada user e joga valor de capacidade para capacidade_real(final)
                idx_mask = find(mask2(user,:));
                real_capacity(user,idx_mask) = capacity(user,idx_mask);
            end
            
        end
   
        %b = sum(capacity);
        b = sum(real_capacity);

        
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
     % loading dynamic max vazao
     DynamicData=load('dynamicMaxVazao_num3.mat');
     D3 = DynamicData.Dynamic.DataSNR;
     D4 = DynamicData.Dynamic.DataBPRB;
end  

%% Gera graficos de Bits/SNR
figure;
plot(SNR, bits_per_rb, '-ok','LineWidth',1.5);
%title('Alocação de Recursos em sistema de multiplo acesso Ortogonal');
xlabel('SNR [dB]'); 
ylabel('Bits/RB'); 
grid on;
grid minor;

hold on;
plot(D1, D2, '--r');
hold on;
plot(D3, D4, '--b');
legend('Alocação Dinâmica','Alocação Estática', 'Dinamica - Máx. Vazão')