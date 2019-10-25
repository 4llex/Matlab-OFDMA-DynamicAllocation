function [ Hrb ] = rb_h_media( H, SC_per_RB )
% INPUTS: H -> Frequency Responcy for each subcarrier
%         X -> Subcarriers per Resource Block
%
% OUTPUT: Hrb -> A vector with the Frequency Response's avarages 
% for 132 Resource Blocks.
%
% Calculate the avarage of an X(depends numerology) amount of subcarriers
% to represents a Resource Block 
    
    qtd_RB = 132;     % Number of Resource Blocks
    RB = 1;
    Hrb = zeros(1,qtd_RB);

    for i=1:SC_per_RB:qtd_RB*SC_per_RB

        Hrb(RB) = sum(H(i:RB*SC_per_RB))/SC_per_RB; %abs??

        RB=RB+1;

    end


end

