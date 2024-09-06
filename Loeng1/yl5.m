clc
clear

F = @(x,y) sin(x-y)/x^2.0 + ( cos(2.0*x+y)/(x-y)^4.0 )^(1.0/3.0);

F(50, -30)

%{
ans =

    0.0021
%}