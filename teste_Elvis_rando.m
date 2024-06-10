%% Algoritmo de lagrangeano aumentado para o problema da mochila 
%% Problema teste 2 Variaveis Randomly generated tests??

clear all;
clc
n=10000000; % ate 10k
if n>1e7
    disp('numero de variaveis ultrapassa a memoria do pc')
    break
end
rep=round((1e7)/n)
pause
iter=zeros(1,rep);
tiempo=zeros(1,rep);
for kk=1:rep
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prec = 1e-4;
%p = ones(n,1);
p=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];
%b = p;
b=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];
%a = (2:n+1)';
a=10 + 15.*rand(n,1); % p=10 + (25-10).*rand(10,1) U[10,25];

B=1 + 14.*rand(2,n);  % B=1 + 14.*rand(2,n) U[1,15];

l=min(B)';
u=max(B)';
btl=l'*b;
btu=u'*b;
c=btl+(btu-btl)*rand(1);

disp('gerou dados');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); % Ponto inicial
lamb=10;        % Multiplicador inicial a restricao de igualdade b'*x-c
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
r=1/(n*10);          % Parametro de penalidade inicial
%r=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=cputime;     % Time
for k = 1:600000

alfa=r/(1+r*(b'*(b./p)));
%lar=a-(lamb-r*c)*b;
%xb=lar./p-(alfa*((b./p)'*lar)*(b./p));

lar=lamb-r*c;

xb=(a-lar*b)./p-r*alfa*(b./p)*((b'*(a./p))-(lar*(b'*(b./p)))); % xb=solucao do grad(L)=0 com r. 
%xb=(a-lar*b)./p-alfa*(b./p)*((b'*(a./p))-(lar*(b'*(b./p))));    % xb=solucao do grad(L)=0.

% gradkk=norm(p.*xb+r*b*(b'*xb)+(lamb-r*c)*b-a)



 %% Proje????o da solucao na caixa
 xk=max(l,min(xb,u));
 
   %% Atualizamos o multiplicador da restricao de igualdade
    rest=b'*xk-c;
    lamb=lamb+r*rest;
      k=k+1;
      
      %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2=abs(rest);
    %pause(0.5)
    if erro1 <= prec && erro2 <= prec
        break
    else
      xk0=xk;
       %if erro2 > 0.5*errot && k1>50
       %      r=10*r
       %     k1=1;
      %else
      %  k1=k1+1;
       %end
      %r=10;
    end
end
timen=cputime-tt; % Tempo utilizado para resolver o problema
k;                % Numero de iteracoes utilizadas para resolver o problema
iter(1,kk)=k;
tiempo(1,kk)=timen;
%pause
end

maxiterr=max(iter)
miniterr=min(iter)
mediaiter=mean(iter)

maxtiempo=max(tiempo)
mintiempo=min(tiempo)
mediatiempo=mean(tiempo)

%% teste da condicao de parada
%gradl=diag(p)*xk+r*b*(b'*xk)+(lamb-r*c)*b-a;
%aux=xk-gradl;
%proj=max(l,min(aux,u));
%test=norm(proj-xk)


%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);