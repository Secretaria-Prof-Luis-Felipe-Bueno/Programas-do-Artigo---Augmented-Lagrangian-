clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resolve o problema [x y z][1 -1 0; -1 2 0; 0 0 1][x y z] sa 
%z=0, [x y z] in [1 0 -1, 5 1 1]

prec = 1e-4;
%p = ones(n,1);

%% gracao de dados
% p= matriz do problema]
n=20;

l = -(n*ones(n,1))/2;
u = (n*ones(n,1))/2;
c = n; b = ones(n,1); 

%gerando P posa defa aleatoria centrada no zero e somada com um multiplo
%gde da identidade, p ficar longe do zero
grande = 1;
numeromenorquen = 20;
pqn = 0.5-rand(numeromenorquen,numeromenorquen);
pqn = pqn'*pqn;
P = grande*eye(n,n);
P(1:numeromenorquen,1:numeromenorquen) = pqn;
a = 0.5-rand(n,1);

%b = p;

% p=10 + (25-10).*rand(10,1) U[10,25];
%a = (2:n+1)';

% p=10 + (25-10).*rand(10,1) U[10,25];

btl=l'*b;
btu=u'*b;

disp('gerou dados')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); %1-rand(n,1); % Ponto inicial
lamb=10;        % Multiplicador inicial a restricao de igualdade b'*x-c
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
r=1;          % Parametro de penalidade inicial
%r=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=cputime; k=1;  aux=0;  % Time
%calculando P inversa b:
pinvb = P\b; pinva = P\a; rest0 = b'*xk0 - c; rest=rest0;
while (k<6000)&&(aux==0)
  
% %   bt=sqrt(r)*b;
% %   alfa=1/(1+bt'*bt);
% %   IP= eye(n)-(alfa*(bt*bt')); % inversa de P formula de Sherman-Morrisson
% %   xk=IP*(a-((lamb-r*c)*b)); % Solucao de grad(lag)=0

alfa=r/(1+r*b'*pinvb);
lar = r*c - lamb;
xk = pinva + (lar - alfa*(b'*pinva + lar*(b'*pinvb)))*pinvb;

xk=max(l,min(xk,u));

    rest=b'*xk-c;
    
    
    if abs(rest)>0.1*abs(rest0)
        r=2*r;
    end   
    rest0 = rest;
        
    %r=n/(4*k+1);x0
    
    %lamb0=lamb;
    lamb=lamb+r*rest;
    %r=abs(lamb-lamb0)/abs(rest0)
      k=k+1;
      
      %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2=abs(rest);
    if erro1 <= prec && erro2 <= prec
        aux=1;
    else
      xk0=xk;  
      %r=10;
    end
    k=k+1;
    %pause
end
timen=cputime-tt
lamb
%xk
k
%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);

%comparando quadprog:
xquadprog = quadprog(P,-a,[],[],b',c,l,u);

