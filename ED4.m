%% Incial
% Pr e Ec
Pr_values = logspace(0, 3, 1000); % Pr variando de 1 a 10^3
Ec_values = [-1, 0, 1]; % Valores de Ec

% Vetor para armazenar as soluções
D_solutions = zeros(length(Pr_values), length(Ec_values));

%% Resolver o calculo numerico
% For do Ec
for j = 1:length(Ec_values)
    Ec = Ec_values(j);

    % For do Pr
    for i = 1:length(Pr_values)
        Pr = Pr_values(i);
        
        % Funçao principal
        F = @(D) ((14 - D.^2) .* D.^3 .* Pr) / 13 - ...
                 (1 + ((3/2) * Pr * Ec * D.^2 .* (1 - (3/2) * D.^2 + (1/5) * D.^4)));
        
        % Resolver numericamente para D usando fzero (considerando um chute inicial)
        D_initial_guess = 0.5; % Chute inicial
        D_solutions(i, j) = fzero(F, D_initial_guess);
    end
end

%% Outro
% Calculando a linha adicional D = Pr^(-1/3)
D_additional = Pr_values.^(-1/3);

%% Plotando os resultados
figure;
hold on;

% for para cada Ec
for j = 1:length(Ec_values)
    plot(Pr_values, D_solutions(:, j), 'DisplayName', ['Ec = ', num2str(Ec_values(j))]);
end

plot(Pr_values, D_additional, '--', 'DisplayName', 'D = Pr^{-1/3}', 'LineWidth', 2);

set(gca, 'XScale', 'log'); % Escala logarítmica para o eixo x (Pr)
set(gca, 'YScale', 'log'); % Escala logarítmica para o eixo y (D)
xlabel('Pr');
ylabel('D');
title('Solução numérica de D em função de Pr para diferentes valores de Ec');
legend('Location', 'best');
grid on;
hold off;

