function sol = run_SIRV()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    

    function dydt = sirv(t,y,Z)
        % SIRV funcrion definition 
        
        S = y(1); I = y(2); R = y(3); V = y(4);
        
        % Calcul de sum_{i = 1}^9 (1-m-beta I)^i V(t - i)
        reimmun = 0 ;
        for i  = 1:(t_immun-1)
            if (t < (i+1))
              reimmun = reimmun + (1-m-beta*I)^i * v*m*N ; % V(t-i-1) remplacé par V à l'équilibre
            else
              reimmun = reimmun + (1-m-beta*I)^i * Z(4,i+1) ; % Z(4,i) = V(t-i-1)
            end
        end

        % Calcul de V(t - t_immun)
        if t < t_immun
            Vlag = v*m*N ; % condition initale pour V
        else
            Vlag = Z(4,t_immun) ; % valeur de V(t - t_immun)
        end

        % Équations
        dydt = zeros(4,1);
        dydt(1) = (1-v)*m*N - m*S - beta*I*S + (1-m-beta*I)^t_immun * Vlag ; % equation de S
        dydt(2) = beta*S*I - m*I - (g-m)*I ; % equation de I
        dydt(3) = (g-m)*I - m*R ; % equation de R
        dydt(4) = v*m*N - m*V - (1-m-beta*I)*V + beta*I*reimmun ; % equation de V
      
    end % end of nested function sirv

% Paramètres du modèle
N = 10 ;      % popultation totale
m = 1/80 ;    % taux de mortalité/natalité, essperance de vie de 80ans
v = 0.8 ;    % couverture vaccinale de 80%RUN_SEIR simulation of the SEIR model
g = 52/3 ;    % durée de l'infection 3semaines
R0 = 6.5 ;      % taux de reproduction de base
beta = R0*(m+g)/N ;   % taux d'infection S -> I
t_immun = 100 ;    % durée de l'immunité vaccinale

  
% Paramètres d'intégration
Si = (1-v)*N + (1-m)^10 * v*N -1 ;
Ii = 1 ;
Ri = 0 ;
Vi = v*m*N;
IC = [Si ; Ii ; Ri ; Vi]; % conditions initales
tspan = [0:150]; % en années
lags = [1:t_immun];

% simulations
%opts = ddeset('RelTol',1e-10,'AbsTol',1e-10);
sol = dde23(@sirv,lags,IC,tspan);

% Affichage 
f1 = figure(1); clf;
plot(sol.x, sol.y(1,:));
title('Évolution de S avec le modèle SIRV');
xlabel('time t');
ylabel('S(t)');

f2 = figure(2); clf;
plot(sol.x, sol.y(2,:));
title('Évolution de I avec le modèle SIRV');
xlabel('time t');
ylabel('I(t)');

f3 = figure(3); clf;
plot(sol.x, sol.y(3,:));
title('Évolution de R avec le modèle SIRV');
xlabel('time t');
ylabel('R(t)');

f4 = figure(4); clf;
plot(sol.x, sol.y(4,:));
title('Évolution de V avec le modèle SIRV');
xlabel('time t');
ylabel('V(t)');
    

    

end

