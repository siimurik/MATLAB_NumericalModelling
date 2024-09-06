clc
clear

A = [3.0 12.0 52.0
    4.0 6.0 -11.0
    -2.0 7.0 2.0 ]

b = [13 -2 5]

display("Ülesanne a)")
B = 2 .* transpose(A) + inv(A)
C = A*B

display("Ülesanne b)")
det(A)

display("Ülesanne c)")
b'

display("Ülesanne d)")
D = A.*B

%{
    "Ülesanne a)"


B =

    6.0354    8.1352   -4.1765
   24.0056   12.0437   14.0958
  104.0159  -22.0179    3.9881


C =

   1.0e+03 *

    5.7150   -0.9760    0.3640
   -0.9760    0.3470    0.0240
    0.3640    0.0240    0.1150

    "Ülesanne b)"


ans =

        2515

    "Ülesanne c)"


ans =

    13
    -2
     5

    "Ülesanne d)"


D =

   18.1062   97.6223 -217.1801
   96.0223   72.2624 -155.0541
 -208.0318 -154.1252    7.9761
%}