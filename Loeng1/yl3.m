clc
clear

% Funktsiooni defineerimine
y=@(x) (x^2.0-3.0)*(2.0+x)^4.0 - 5.0*exp(x) + 2*cos(x+1.0);

y(-1)
y(4)
y(3)

%{
ans =

   -1.8394


ans =

   1.6576e+04


ans =

   3.6483e+03
%}