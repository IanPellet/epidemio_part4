function sol = run_SIRV()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    function dydt = sirv(t,y,Z)
        % SIRV funcrion definition 
        
        S = y(1); I = y(2); R = y(3); V = 0;
        
        % Calcul de V(t - t_immun)
        %if t < t_immun
        %    Vlag = v*m*N ; % condition initale pour V
        %else
        %    Vlag = Z(4,t_immun) ; % valeur de V(t - t_immun)
        %end
        
        %somme_dV1 = 0 ;
        %for i = 0:t_immun
        %    if (t < (i+1))
        %        somme_dV1 = somme_dV1 + (1-m-beta*I)^i * v*m*N ; % V(t-i-1) remplacé par V à l'équilibre
        %    else
        %        somme_dV1 = somme_dV1 + (1-m-beta)^i * Z(4,i+1) ; % Z(4,i) = V(t-i-1)
        %    end
        %end        
        Vlag = 0;
        % Équations
        dydt = zeros(3,1);
        %dydt(1) = (1-v)*m*N - m*S - beta*I*S + (1-m-beta*I)^t_immun * Vlag ; % equation de S
        dydt(1) = (1-v)*m*N - m*S - beta*I*S  ; % equation de S
        dydt(2) = beta*S*I - m*I - g*I ; % equation de I
        dydt(3) = g*I + v*m*N - m*R ; % equation de R
        %dydt(4) = v*m*N + beta*I*somme_dV1 - V ; % equation de V
        
        %disp(dydt)
      
    end % end of nested function sirv

% Paramètres du modèle
N = 1e06 ;      % popultation totale
m = 1/80 ;    % taux de mortalité/natalité, essperance de vie de 80ans
v = 0.0 ;    % couverture vaccinale de 80%RUN_SEIR simulation of the SEIR model
g = 52/3 ;    % durée de l'infection 3semaines
R0 = 6.5 ;      % taux de reproduction de base
beta = R0*(m+g)/N ;   % taux d'infection S -> I
t_immun = 1000 ;    % durée de l'immunité vaccinale



  
% Paramètres d'intégration
%Si = (1-v)*N + (1-m)^10 * v*N  ;
Si = (1-v)*N -1 ;
Ii = 1 ;
Ri = v*N ;
%Vi = v*m*N;
IC = [Si ; Ii ; Ri ]; % conditions initales
tspan = [0:0.001:200]; % en années
lags = [1:20];
% simulations
%opts = ddeset('RelTol',1e-10,'AbsTol',1e-10);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = dde23(@sirv,lags,IC,tspan,options);

% Affichage 
f1 = figure(4); clf;
plot(sol.x, sol.y(1,:));
title('Évolution de S avec le modèle SIRV');
xlabel('time t');
ylabel('S(t)');

f2 = figure(5); clf;
plot(sol.x, sol.y(2,:));
title('Évolution de I avec le modèle SIRV');
xlabel('time t');
ylabel('I(t)');

f3 = figure(6); clf;
plot(sol.x, sol.y(3,:));
title('Évolution de R avec le modèle SIRV');
xlabel('time t');
ylabel('R(t)');



    

end

