function sol = SIR()
% Modèle SIR avec vaccination

% Paramètres du modèle
N = 1e6 ;      % popultation totale
m = 1/80 ;    % taux de mortalité/natalité, essperance de vie de 80ans
v = 0.8 ;    % couverture vaccinale de 80%RUN_SEIR simulation of the SEIR model
g = 52/3 ;    % durée de l'infection 3semaines
R0 = 6.5 ;      % taux de reproduction de base
beta = R0*(m+g)/N ;   % taux d'infection S -> I
% t_immun = 10 ;    % durée de l'immunité vaccinale

  
% Paramètres d'intégration
Si = (1-v)*N - 1 ;
Ii = 1 ;
Ri = v*N ;
% Vi = v*m*N;
% times = seq(0, 2, 0.1);
IC = [Si ; Ii ; Ri]; % seed a few exposed individuals
tspan = [0:0.001:200]; % en années

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options);

%figure(2); clf;
% semilogy(sol.x, sol.y(1:3,:));
f1 = figure(1); clf;
plot(sol.x, sol.y(1,:));
f2 = figure(2); clf;
plot(sol.x, sol.y(2,:));
f3 = figure(3); clf;
plot(sol.x, sol.y(3,:));

% disp(['Percentage of contaminated 100*(E+I+R)/N at day ' num2str(tspan(2)), ': ', ...
 % num2str(100*sum(sol.y(2:3,end)/N))]);


    function dxdt = sir(~,xx)  % nested function 
        % SIR is the ODE rhs 
        
        S = xx(1); I = xx(2); R = xx(3);
        
        dxdt = [ (1-v)*m*N - m*S - beta*S*I ;
                 beta*S*I - m*I - g*I ;
                 g*I + v*m*N - m*R ]; 
    end                         % end nested function seir

end                             % end main function run_SEIR 
