clear;
clc;
N = 3;
%
[K,g] = gera_LCP(N);
%
lamb = 1; % multiplicador de lagrange inicial
r = 1000; % parametro de penalidade (fixo)
gama=1;     % Parametro de regularizacao do lagrangeano
xk = ones(N,1); % ponto inicial
acur = 1e-4; % acuracia
kmax = 6000;
% pkg load optim
pq  = quadprog(2*K,g,-K,g,[],[],zeros(N,1),[])
pqt = quadprog(K,g,g',0,[],[],zeros(N,1),[]) 
[kgd, jgd, tgd, e1gd,e2gd,e3gd,sol,s] = Mochila_lcp_1(xk,K,g,N,kmax,acur,lamb,r,gama,'normal');
%[kgd, jgd, tgd, e1gd,e2gd,e3gd] = gradiente_lcp_1(xk,K,g,N,kmax,acur,lamb,r);
fk2 = abs(sol'*(K*sol + g))
e22 = g'*sol + s
fk = abs(pq'*(K*pq + g))
e1 = norm(min(pq,0))
e2 = norm(min((K*pq + g),0))
sol