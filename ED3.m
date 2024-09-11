clear; clc;

ED()

function ED()
    %% Parâmetros do problema

    n = 1000;  % Nós
    mu = 0.03; % Viscosidade (kg/ms)
    k = 0.6;   % Condutividade (W/mK)

    l = 0.1;   % Altura (m)
    Tw = 298;  % Temperatura ambiente em Kelvin
    U = 10; % Velocidade (m/s)
    G = -0.5;  % Dp/dx

    
    %% Resolução
    u = sol_u(l, U, G, mu, n);
    T = sol_T(l, Tw, mu, k, u, n);
    
    % Eixo y 
    y = linspace(0, l, n);
    
    %% Soluções analíticas
    % Solução analítica para velocidade e temperatura
    eta = y / l;
    Up = (l^2 / (12 * mu)) * (-G);
    
    % Chama as soluções analíticas para armazená-las
    u_analitica = sol_analitica_u(eta, U, Up);
    T_analitica = sol_analitica_T(eta, U, Up, Tw, mu, k);
    
    %% Plotagem
    
    % Gráfico de velocidade
    figure;
    plot(u, y, 'r', 'DisplayName', 'Numerico');
    hold on;
    plot(u_analitica, y, '--b', 'DisplayName', 'Analitico');
    xlabel('Velocidade [m/s]');
    ylabel('Altura [m]');
    title('Velocidade');
    legend;
    grid on;
    
    % Gráfico de temperatura
    figure;
    plot(T, y, 'r', 'DisplayName', 'Numerico');
    hold on;
    plot(T_analitica, y, '--b', 'DisplayName', 'Analitico');
    xlabel('Temperatura [K]');
    ylabel('Altura [m]');
    title('Temperatura');
    legend;
    grid on;
end


%% Funções para resolução

function x = TDMA(a, b, c, d)
    n = length(d);
    % Eliminação para frente
    for i = 2:n
        m = a(i-1) / b(i-1);
        b(i) = b(i) - m * c(i-1);
        d(i) = d(i) - m * d(i-1);
    end
    % Substituição para trás
    x = zeros(n, 1);
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end

function u = sol_u(l, U, G, mu, n)
    h = l / (n + 1);
    % Coeficientes do sistema tridiagonal
    a = ones(n-1, 1);
    b = -2 * ones(n, 1);
    c = ones(n-1, 1);
    d = (G * h^2 / mu) * ones(n, 1);
    
    % Condições de contorno
    b(1) = 1; c(1) = 0; d(1) = 0;
    a(end) = 0; b(end) = 1; d(end) = U;
    
    % Solução de u
    u = TDMA(a, b, c, d);
end

function T = sol_T(l, Tw, mu, k, u, n)
    h = l / (n + 1);
    % Coeficientes do sistema tridiagonal
    a = ones(n-1, 1);
    b = -2 * ones(n, 1);
    c = ones(n-1, 1);
    d = zeros(n, 1);
    
    % Condições de contorno T(0) = Tw
    b(1) = 1; c(1) = 0; d(1) = Tw;
    
    % Equações de diferença para T
    for i = 2:n-1
        d(i) = -mu / k * (u(i+1) - u(i))^2; 
    end
    
    % Condição de contorno para T'(l) = 0
    b(n) = b(n) + c(n-1);  
    d(n) = -mu / k * (u(n) - u(n-1))^2; 
    
    % Solução de T
    T = TDMA(a, b, c, d);
end

function u_analitica = sol_analitica_u(eta, U, Up)
    u_analitica = U * eta + 6 * Up * eta .* (1 - eta);
end

function T_analitica = sol_analitica_T(eta, U, Up, Tw, mu, k)
    T_analitica = ((mu * U^2) / (2 * k)) * (eta .* (2 - eta) ...
        - 4 * Up / U * eta.^2 .* (3 - 2 * eta) ...
        + 12 * (Up / U)^2 * eta .* (2 - 3 * eta + 4 * eta.^2 - 2 * eta.^3)) + Tw;
end

