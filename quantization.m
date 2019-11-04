function [capacity_quantizada] = quantizar(b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    b(b<2) = 0;
    b((b>2)&(b<4)) = 2;
    b((b>4)&(b<6)) = 4;
    b((b>6)&(b<8)) = 6;
    b(b>8) = 8;

    capacity_quantizada = b;
end

