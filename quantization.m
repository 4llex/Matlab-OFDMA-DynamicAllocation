function [capacity_quantizada] = quantization(b, Num)
%   UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    b(b<0,083) = 0;
    b((b>=0.083)&(b<0.167)) = 0.083;
    b((b>=0.167)&(b<0.250)) = 0.167;
    b((b>=0.250)&(b<0.333)) = 0.250;
    b((b>=0.333)&(b<0.417)) = 0.333;
    b((b>=0.417)&(b<0.583)) = 0.417;
    b((b>=0.583)&(b<0.750)) = 0.583;
    b((b>=0.750)&(b<0.833)) = 0.750;
    b((b>=0.833)&(b<1.000)) = 0.833;
    b((b>=1.000)&(b<1.166)) = 1.000;
    b((b>=1.166)&(b<1.500)) = 1.166;
    b((b>=1.500)&(b<1.833)) = 1.500;
    b((b>=1.833)&(b<2.167)) = 1.833;
    b((b>=2.167)&(b<2.500)) = 2.167;
    b((b>=2.500)&(b<3.000)) = 2.500;
    b((b>=3.000)&(b<3.333)) = 3.000;
    b((b>=3.333)&(b<3.500)) = 3.333;
    b((b>=3.500)&(b<4.000)) = 3.500;
    b((b>=4.000)&(b<4.500)) = 4.000;
    b((b>=4.500)&(b<4.750)) = 4.500;
    b((b>=4.750)&(b<5.250)) = 4.750;
    b((b>=5.250)&(b<5.500)) = 5.250;
    if (Num == 3)
        b(b>=5.500) = 5.500;
    else 
        b((b>=5.500)&(b<6.000)) = 5.500;
        b((b>=6.000)&(b<6.667)) = 6.000;
        b((b>=6.667)&(b<7.000)) = 6.667;
        b((b>=7.000)&(b<7.333)) = 7.000;
        b((b>=7.333)&(b<7.667)) = 7.333;
        b(b>=7.667) = 7.667;    
    end

    capacity_quantizada = b;
end

