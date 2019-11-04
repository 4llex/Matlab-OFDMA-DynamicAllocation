function [ subPower, subAlloc, subCapacity ] = fcn_waterfilling(Power,SigmaSqr,Gamma,H,subAlloc)
% -------------------------------------------------------------------------
% This function distributes the specified amount of power among the
% subcarriers using the waterfilling algorithm
% -------------------------------------------------------------------------
%   Power    : power to be distributed among the subcarriers
%   SigmaSqr : noise variance
%   subGains : vector with tha channel gain at each subcarrier
%   TargetSer: target ser
%   subAlloc : boolean vector indicating which subcarriers are allocated
%   subPower : vector indicating the power allocated to each subcarrier



%AllowedSubCarriers=subAlloc;
subPower=ones(size(subAlloc))*-1;
while(sum(subPower<0))
    subPower(subPower<0)=0;
    idx = find(subAlloc>0);
    subPower(idx)=(Power/sum(subAlloc))...
        + ((Gamma*SigmaSqr)/sum(subAlloc))*sum(abs(H(idx)).^-2)...
        - (Gamma*SigmaSqr)*(abs(H(idx)).^-2);
    subAlloc(subPower<0)=0;
end
%subAlloc=AllowedSubCarriers;
%subCapacity = log2(1 + (subPower.*(abs(H).^2))/(Gamma*SigmaSqr));
subCapacity = log2(1 + (subPower.*(abs(H).^2))/(Gamma*SigmaSqr));  % Modulation constraint
