0;
% ========================================
% Authors: 
% - 9298614 - Arthur Prado de Fazio
% - 11796378 - João Henri Carrenho Rocha
% ========================================

function print_iteration_info(c, bind, x, iteration)
    printf("Iterando %d\n", iteration);
    for idx = 1:numel(bind)
        printf("%d %f\n", bind(idx), x(bind(idx)));
    endfor

    printf("\nValor função objetivo: %f\n", c'*x);

endfunction

function [reduced_costs, non_basic_variables] = calculate_reduced_costs(A, c, n, bind, Binv)
    non_basic_variables = setdiff(1:n , bind);
    reduced_costs = c(non_basic_variables)' - c(bind)' * Binv * A(:,non_basic_variables);
    
    printf("\nCustos reduzidos\n");
    for idx = 1:numel(reduced_costs)
        printf("%d %f\n", non_basic_variables(idx), reduced_costs(idx) );
    endfor
endfunction

function u = compute_u(A, bind, Binv, j)
    u = Binv * A(:, j);
    
    printf("\nDireção\n");
    for idx =  1:numel(u)
        printf("%d %f\n", bind(idx), u(idx))
    endfor
endfunction

function [theta, l] = calculate_theta(m, x, bind, u)
    theta = Inf;
    l = -1;
    
    for i = 1:m
        if (u(i) > 0)
            current_theta = x(bind(i))/u(i);
            if (current_theta < theta)
                l = i;
                theta = current_theta;
            endif
        endif
    endfor
endfunction

function [ind v] = optimal_solution(c, n, x, bind)
    ind = 0;
    v = x;

    printf("\nSolução ótima encontrada com custo %f\n", c'*x);
    for i = 1:numel(v)
        printf("%d %f\n", i, v(i));
    endfor
endfunction

function [ind v] = unbounded_direction(u, bind, j, m, n)
    ind = -1;
    v = zeros(n, 1);
    v(j) = 1;
    for i = 1:m
      v(bind(i)) = -u(i);
    endfor

    printf("\nProblema é ilimitado. Direção de custo -Inf:\n");
    for i = 1:numel(v)
        printf("%d %f\n", i, v(i));
    endfor
  endfunction

function x = update_basic_solution(x, bind, m, u, theta, j, l)
    x(j) = theta;
    for i = 1:m
        if (i != l)
            x(bind(i)) -= theta*u(i);
        else
            x(bind(i)) = 0;
        endif
    endfor
 endfunction

 function Binv = update_Binv(Binv, m, u, l)
    for i = 1:m
        if (i != l)
            Binv(i, :) += Binv(l, :)*(-u(i)/u(l));
        endif
    endfor

    Binv(l, :) = Binv(l, :) / (u(l));
 endfunction

function [ind v ] = simplex(A,b,c,m,n,x,bind, Binv)
    iteration = 0;

    while(true)
        print_iteration_info(c, bind, x, iteration)

        [reduced_costs, non_basic_variables] = calculate_reduced_costs(A, c, n, bind, Binv);

        if (all(reduced_costs >= 0))   
            [ind v] = optimal_solution(c, n, x, bind);     
            return
        endif

        j = non_basic_variables(find(reduced_costs < 0)(1));
        printf("\nEntra na base: %d\n", j)

        u = compute_u(A, bind, Binv, j);
        [theta, l] = calculate_theta(m, x, bind, u);        
        printf("\nTheta*\n%f\n", theta);
        
        if (theta == Inf)
            [ind v] = unbounded_direction(u, bind, j, m, n);
            return
        endif

        printf("\nSai da base: %d\n\n", bind(l));

        x = update_basic_solution(x, bind, m, u, theta, j, l);
        bind(l) = j;
        Binv = update_Binv(Binv, m, u, l);
        iteration++;
    endwhile
endfunction

% Example with optimal solution:
% b = [20 20 20]';
% c = [-10  -12 -12 0 0 0]';
% A = [ [1; 2; 2] [2; 1; 2] [2; 2; 1] [1; 0; 0] [0; 1; 0] [0; 0; 1] ];
% m = 3;
% n = 6;
% x = [0, 0, 0, 20, 20, 20]';
% bind = [4, 5, 6];
% B = [ [1; 0; 0] [0; 1; 0] [0; 0; 1]];
% Binv = inv(B);

% Example with unlimited cost:
A = [ [1; 0] [2; 1] [0; 1] [1; 1] [0; 0] ];
b = [10 3]';
c = [4  5  1 -1 -1]';
m = 2;
n = 5;
x = [10 0 3 0 0]';
bind = [1 3];
B = [ [1; 0] [0; 1] ];
Binv = inv(B);

[ind v] = simplex(A,b,c,m,n,x,bind, Binv);
