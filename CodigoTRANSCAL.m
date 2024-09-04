%% PlacaVeleTemp 

% FALTANDO b_T'

%% VELOCIDADE _______________________________________________________________________
n = 10;         % Malha
l = 1;           % Comprimento 
h = l/n;         % Passo
mu = 100;  % Viscosidade 
k = 1000;        % Constante
U = 5;          % Velocidade final
T_w = 100;       % Temperatura na parede

%vetor de velocidade
u = zeros(1, n);

% Construção da matriz A e b
Au = zeros(n);
bu = zeros(n, 1);

%%
for i = 1:n
    if i == 1  % Condição de Contorno Inicio
        Au(i, i) = 1;
        Au(i, i+1) = 0;
        bu(i) = 0;
    elseif i == n  % Condição de Contorno Fim
        Au(i, i-1) = 0;
        Au(i, i) = 1;
        bu(i) = U;
    else
        Au(i, i-1) = 1;
        Au(i, i) = -2;
        Au(i, i+1) = 1;
        bu(i) = (k*(h^2))/mu;
    end
end
%%

u = inv(Au)*bu;


% % Eliminação de Gauss para resolver o sistema tridiagonal
% for i = 2:n
%     factor = Au(i, i-1) / Au(i-1, i-1);
%     Au(i, :) = Au(i, :) - factor * Au(i-1, :);
%     bu(i) = bu(i) - factor * bu(i-1);
% end
% 
% % Substituição inversa
% u(n) = bu(n) / Au(n, n);
% for i = n-1:-1:1
%     u(i) = (bu(i) - Au(i, i+1) * u(i+1)) / Au(i, i);
% end

%% Temperatura ____________________________________________________________________________
% Inicialização do vetor de temperatura
T = zeros(1, n);

% Matriz A e b
A_T = zeros(n);
b_T = zeros(n, 1);
for i = 1:n
    if i == 1
        A_T(i, i) = 1;
        A_T(i, i+1) = 0;
        b_T(i) = T_w;
    elseif i == n
           A_T = 1;
           A_T(i, i-1) = 0;
%         
%         A_T(i, i-1) = 1;
%         A_T(i, i) = -2;
    else
        A_T(i, i-1) = 1;
        A_T(i, i) = -2;
        A_T(i, i+1) = 1;
    end
end


for i = 2:n-1
    b_T(i) = -(mu/(k*4)) * ((u(i+1) - u(i-1))^2);
end


%%

% T = inv(A_T)*b_T;
% 
% % % Resolução da equação
% % for i = 2:n-1
% %     b_T(i) = -(mu/(k*4)) * ((u(i+1) - u(i-1))^2);
% % end
% % 
% % 
% % % Eliminação de Gauss para resolver o sistema tridiagonal
% % for i = 2:n
% %     factor = A_T(i, i-1) / A_T(i-1, i-1);
% %     A_T(i, :) =A_T(i, :) - factor * A_T(i-1, :);
% % end
% % 
% 
% 
% % Substituição inversa
% T(n) = b_T(n) / A_T(n, n);
% for i = n-1:-1:1
%     T(i) = (b_T(i) - A_T(i, i+1) * T(i+1)) / A_T(i, i);
% end
% 
% %% Solução Analitica 
% 
% y = linspace(0, l, n);
% U_p = ((l^2)/(12*mu))*(-k);
% N = y/l;
% 
% u_analitico = U.*(y/l) + 6*U_p.*(y/l).*(1-(y/l));
% T_analitico =( N.*(2-N) - ((4.*U_p)/U).*(N.^2).*(3-2.*N)  +  12.*((U_p/U).^2).*N.*(2-3.*N +4.*(N.^2)- 2.*(N.^3))  ).*((mu.*U.^2)/(2*k)) + T_w;
% 
% 
% 
% 
% 
% 
% %% Plotagem da solução _______________________________________________________________________
% 
% % 
% % plot(y, u, 'b', 'LineWidth', 2); hold on;
% % plot(y, u_analitico, 'r--', 'LineWidth', 2);
% % xlabel('y');
% % ylabel('u');
% % legend('u', 'u_{analitico}');
% % title('Plot de u e u_{analitico} em função de y');
% % grid on;
% 
% 
% 
% 
% plot(y, T, 'b', 'LineWidth', 2); hold on;
% plot(y, T_analitico, 'r--', 'LineWidth', 2);
% xlabel('y');
% ylabel('T');
% legend('T', 'T_{analitico}');
% title('Plot de u e u_{analitico} em função de y');
% grid on;

