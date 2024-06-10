clear; clc;
tipo = 4;
n = 50000;
% carregando o problema 
www=fullfile(['prob_4_50000_1.txt']);
M=(dlmread(www))';
% carregando o problema 
b = M(1,1:n)';                % left side knapsack constraint  
g = M(1,(n+2):(2*n+1))';      % \gamma_j \in U(1,5)
l = M(1,(2*n+2):(3*n+1))';    %  l left side box constraint    
u = M(1,(3*n+2):(4*n+1))';    %  u right side box constraint
c = M(1,(n+1));               % \hat{b}.
% Lagrangeano aumentado
r = 0.8;
lamb = 10;
[xla,tla,kla,ikla] = la_fuel(g,c,l,u,r,lamb);
% 
% 
vla = fun1_fuel(xla,g);
ilafeas = feas_fuel(xla,c,l,u);
ilaopt = solut_fuel(xla,g,c,l,u);
% 
% 
% Bisection method
[xbs,tbs,kbs,A,B] = bi_section_fuel_v1(g,c,l,u);
xbs = xbs';
ibsfeas = feas_fuel(xbs,c,l,u);
ibsopt = solut_fuel(xbs,g,c,l,u);
vbs = fun1_fuel(xbs,g);

% 
% 
% regula falsi method 
[xrf,trf,krf] = regula_falsi_fuel_v1(g,c,l,u);
xrf = xrf';
irffeas = feas_fuel(xrf,c,l,u);
irfopt = solut_fuel(xrf,g,c,l,u);
vrf = fun1_fuel(xrf,g);
% 
% secant method
[xsc,tsc,ksc] = secant_fuel_v1(g,c,l,u);
xsc = xsc';
iscfeas = feas_fuel(xsc,c,l,u);
iscopt = solut_fuel(xsc,g,c,l,u);
vsc = fun1_fuel(xsc,g);

clc;
disp('A comparison with multiplier methods');
disp('Alg   It.   Time       Feas.         Opt.      F(xk) ');
disp('---------------------------------------------------------');
disp(['RF',+'  | ', +num2str(krf),+'  | ',+num2str(trf),+' | ',+num2str(irffeas),+' | ',+num2str(irfopt),+' | ',+num2str(vrf)]);
disp(['SC',+'  | ', +num2str(ksc),+' | ',+num2str(tsc),+'  | ',+num2str(iscfeas),+'  | ',+num2str(iscopt),+' | ',+num2str(vsc)]);
disp(['BS ',+' | ', +num2str(kbs),+' | ',+num2str(tbs),+' | ',+num2str(ibsfeas),+' | ',+num2str(ibsopt),+' | ',+num2str(vbs)]);
disp(['AL ',+' | ', +num2str(kla),+' | ',+num2str(tla),+'  | ',+num2str(ilafeas),+' | ',+num2str(ilaopt),+' | ',+num2str(vla)]);
disp('RF = Regula falsi');
disp('SC = Secant');
disp('BS = Bisection');
disp('AL = Augmented lagrangian');

