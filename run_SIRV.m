function sol = run_SIRV()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    function dydt = sirv(~,y)
        % SIRV funcrion definition 
        
        S = y(1); I = y(2); R = y(3); V0 = y(4);

        % Équations
        dydt = zeros(3+t_immun,1);
        
        % equation de S
        dydt(1) = (1-v)*m*N - m*S - beta*I*S + (1-m-beta*I)*y(4+t_immun-1) ; % y(4+t_immun-1) = V_{t_immun - 1}
        dydt(2) = beta*S*I - m*I - g*I ; % equation de I
        dydt(3) = g*I - m*R ; % equation de R
        % equation de V
        dydt(4) = v*m*N + beta*I*sum(y(4:(4+t_immun-1))) - V0 ; % sum(y(4:(t_immun-1))) = sum(Vi)
        % équation de Vi
        for i = 1:(t_immun-1) 
            dydt(4+i) = (1-beta*I-m) * y(4 + (i-1)) - y(4 + i) ; % y(4 + (i-1)) = V_{i-1}
        end 
        
        %disp(dydt)
      
    end % end of nested function sirv

% Paramètres du modèle
N = 1e06 ;      % popultation totale
m = 1/80 ;    % taux de mortalité/natalité, essperance de vie de 80ans
v = 0.8 ;    % couverture vaccinale de 80%RUN_SEIR simulation of the SEIR model
g = 52/3 ;    % durée de l'infection 3semaines
R0 = 6.5 ;      % taux de reproduction de base
beta = R0*(m+g)/N ;   % taux d'infection S -> I
t_immun = 10 ;    % durée de l'immunité vaccinale



  
% Paramètres d'intégration
%Si = (1-v)*N + (1-m)^10 * v*N  ;
Si = (1-v)*N-1  ;
Ii = 1 ;
Ri = 0 ;
Vi = v*m*N;
% conditions initales
IC = ones(t_immun+3, 1)*Vi ;
IC(1) = Si ;
IC(2) = Ii ;
IC(3) = Ri ;
 
tspan = [0,5000]; % en années
%lags = [1:200-t_immun];
% simulations
%opts = ddeset('RelTol',1e-10,'AbsTol',1e-10);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
%sol = dde23(@sirv,lags,IC,tspan,options);
sol = ode45(@sirv,tspan,IC,options);
x = linspace(0,5000);
y = deval(sol, x);

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

figure(5); clf;
plot(sol.x, sum(sol.y(1:(t_immun+3),:)));
title('Évolution de N avec le modèle SIRV');
xlabel('time t');
ylabel('N(t)');

figure(6); clf;
plot(x,y(5:13,:));
title('Évolution de Vi avec le modèle SIRV');
xlabel('time t');
ylabel('dVi(t)');
legend;
    

end

