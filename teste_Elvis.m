%% Algoritmo de lagrangeano aumentado para o problema da mochila 
%% Problema teste 1 P=Matriz identidade
clear all;
clc
n=5000; % numero de vaiaveis
if n>1e7
    disp('numero de variaveis ultrapassa a memoria do pc');
end
%rep=round((1e7)/n);
rep=1;
%pause(0.5)
iter=zeros(1,rep);
tiempo=zeros(1,rep);
for kk=1:rep
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados do problema  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prec = 1e-4;
prec2= 1e-4;
p = ones(n,1);
b = p;
a = (2:n+1)';
l = zeros(n,1);%rest=n/10*(b'*xk-c);
u = 10*p;
c = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); % Ponto inicial

%lamb=-100;    % Multiplicador inicial a restricao de igualdade b'*x-c
lamb=10;       % multiplicador inicial

%r=(n+sqrt(n^2+4*n))/2; % penalty parameter
%r=1/10;   % I a best penalty parameter with r
r=100;         % Is a best penalty parameter without r=1 , 10, 100 are test 
%r=1/(10*n);            % Is a very good penalty parameter for ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tt=cputime; k=1; aux=0;% Time 
while k <= (60000)&&(aux==0)
        
    %% Passo 1 calculo do gradiente
alfa=r/(1+r*(b'*b));
lar=lamb-r*c;
%xb=a-lar*b-r*alfa*b*((b'*a)-(lar*(b'*b))); % xb=solucao do grad(L)=0 com r.
xb=a-lar*b-alfa*b*((b'*a)-(lar*(b'*b)));    % xb=solucao do grad(L)=0 sem r.



%ngradl=norm(xb+r*b*(b'*xb)+lar*b-a)       % Teste da norma do grad(L)=0
 %% Projeçao da solucao na caixa
 xk=max(l,min(xb,u));
 
 %% Passo 2  
 %% Atualizamos o multiplicador da restricao de igualdade
 
    rest=(b'*xk)-c;
    lamb=lamb+r*rest; 
      
  %% Criterio de parada
    erro1=norm(xk-xk0);  % criterio do ponto xk
    erro2=abs(rest)    % criterio das restricoes
    if erro1 <= prec && erro2 <= prec2
        aux=1; 
    else
      xk0=xk;
    end
    
%gradl=xk+r*b*(b'*xk)+(lamb-r*c)'*b-a;
%aux=xk-gradl;
%proj=max(l,min(aux,u));
%projtest=norm(proj-xk)   
k=k+1;
    
end
timen=cputime-tt; % Tempo utilizado para resolver o problema
k;         % Numero de iteracoes utilizadas para resolver o problema
iter(1,kk)=k;
tiempo(1,kk)=timen;
end %for do rep

maxiterr=max(iter)
miniterr=min(iter)
mediaiter=mean(iter)

maxtiempo=max(tiempo)
mintiempo=min(tiempo)
mediatiempo=mean(tiempo)




%% Teste da solucao |proj(xk-grad(Lxk))-xk|
%gradl=xk+r*b*(b'*xk)+(lamb-r*c)'*b-a;
%aux=xk-gradl;
%proj=max(l,min(aux,u));
%projtest=norm(proj-xk)