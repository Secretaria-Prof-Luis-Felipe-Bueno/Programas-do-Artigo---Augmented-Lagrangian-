prec = 1e-7;

n = 3;


%aqui sao so os dados

  l = -(n*ones(n,1))/2;
  u = (n*ones(n,1))/2;

  c = 0; b = ones(n,1); 

%gerando P posa defa aleatoria centrada no zero e somada com um multiplo
%gde da identidade, p ficar longe do zero
grande = 1;
numeromenorquen = n;
pqn = 0.5-rand(numeromenorquen,numeromenorquen);
pqn = pqn'*pqn;

  P = grande*eye(n,n);
  P(1:numeromenorquen,1:numeromenorquen) = pqn;
  P=P+grande*eye(n,n);
  a = 0.5-rand(n,1);
  
 %comeca a contagem de tempo 
   

btl=l'*b;
btu=u'*b;


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
while (k<10000)&&(aux==0)
  

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

%agora resolvendo quadprog com o rapaz como ponto inicial
%xk0 = quadprog(P,-a,[],[],b',c,l,u,xk0);
[xk0,fk0,~,out]=fmincon(@(x) 0.5*x'*(P*x)-a'*x,xk0,[],[],b',[c],l,u);
out
pause

timen=cputime-tt
lamb
%xk
k
%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);

%comparando quadprog:
tic
xquadprog = quadprog(P,-a,[],[],b',c,l,u);
%[xquadprog,fquad]=fmincon(@(x) 0.5*x'*(P*x)-a'*x,ones(n,1),[],[],b',[c],l,u);
timenprog = toc