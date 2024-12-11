clc
clear
fa1 = @(x) x.^2;     % x <= 0
fa2 = @(x) -x.^2;    % 0 < x <= 1
fa3 = @(x) 1 - 2.*x; % x > 1
funcs = {fa1, fa2, fa3};
points = [0, 1];

checkSplineContinuity(funcs, points);


function checkSplineContinuity(funcs, points)
    fprintf("Ülesanne 1.\n");
    fprintf("a)\n");
    
    % Display continuity verification for the function values
    fprintf("\nKontrollime, kas esitab splaini (kui saame 0, siis ei):\n")
    for i = 1:length(points)
        disp(funcs{i}(points(i)) == funcs{i+1}(points(i)));
    end
    fprintf("Kõik on tõesed seega edaspidi leiame, mis on siledusaste.\n")

    syms x;
    
    % Calculate first and second derivatives for each function
    first_derivatives = cell(1, length(funcs));
    second_derivatives = cell(1, length(funcs));
    
    for i = 1:length(funcs)
        % Differentiate the functions symbolically and convert to function handles
        first_derivatives{i} = matlabFunction(diff(funcs{i}(x), x));
        second_derivatives{i} = matlabFunction(diff(first_derivatives{i}(x), x));
    end
    
    % Check first derivative continuity (p = 1)
    fprintf("Kontrollime võrdsust, et leida siledusaste.\nKontroll, kas p=1:\n")
    for i = 1:length(points)
        disp(first_derivatives{i}(points(i)) == first_derivatives{i+1}(points(i)));
    end
    
    % Check second derivative continuity (p = 2)
    fprintf("\nKontroll, kas p=2:\n")
    for i = 1:length(points)
        disp(second_derivatives{i}(points(i)) == second_derivatives{i+1}(points(i)));
    end
    
    fprintf("Pole võrdsed, seega jääb kehtima viimane aste.\n");
end
