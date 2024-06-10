function [k,timen,xk0,lamb,r,epslon,xquadprog,timeprog] = func_teste_nda_aleatorio(flag,P,a,b,c,l,u,lambda,pareps,pardesc)

prec = 1e-4;

%xquadprog = quadprog(P,-a,[],[],b',c,l,u);
%flags:
%0: poe os dados default com p randomica em vez de tomar dos inputs
%1 em diante: pega dados do input pra passar pra funcao main
%2: termina qd r fica mt grande e decide usar o solver instead
%3: qd r eh mmuito grande usa o solver so pra resolver o subproblema

%p = ones(n,1);

%% gracao de dados
% p= matriz do problema]
if flag==0
  n = 200;
else
  n = length(a);  
end

if (nargin<6)||(flag==0)
  l = -(n*ones(n,1))/2;
  u = (n*ones(n,1))/2;
end

if (nargin<4)||(flag==0)
  c = 0; b = ones(n,1); 
end

if nargin<9
    pareps = prec; %comeca aqui
    pardesc = 0; %e decai esse tanto a cada vez
end
if nargin<10
   pardesc = 0; 
end

%%
% flag 0 : problema default
% flag 1 : alg normal
% flag 2 : troca aux para r>10000, aux=2 -> quadprog
% flag 3 : hibrido com reset de r, vai pra flag 5
% flag 4 : hibrido mas incorpora erro3 no criterio de parada
% flag 5 : troca para o LA
% flag 6 : velho, troca o erro epslon a cada iter
% flag 7 : velho, troca o erro epslon a cada iter
% flag 8 : velho, troca o erro epslon a cada iter
% flag 9 : erro4, esquema do lambda velho, com matriz
% flag 10 : erro4, esquema do lambda
% flag 11 : usa outro metodo (regula falsi, etc)


%gerando P posa defa aleatoria centrada no zero e somada com um multiplo
%gde da identidade, p ficar longe do zero
%grande = ones(n,1);
grande = zeros(n,1); for i=1:n grande(i)=i; end

numeromenorquen = 200;
pqn = 0.5-rand(numeromenorquen,numeromenorquen);
pqn = pqn'*pqn;

if flag==0
  P = diag(grande);%grande*eye(n,n);
  P(1:numeromenorquen,1:numeromenorquen) = pqn;
  a = 0.5-rand(n,1);
end

%b = p;

% p=10 + (25-10).*rand(10,1) U[10,25];
%a = (2:n+1)';

% p=10 + (25-10).*rand(10,1) U[10,25];

btl=l'*b;
btu=u'*b;

disp('gerou dados');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dados inciais  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xk0=ones(n,1); %1-rand(n,1); % Ponto inicial
if nargin <8
  lamb = 2; %lamb=10;        % Multiplicador inicial a restricao de igualdade b'*x-c
else
  lamb = lambda;  
end
v=ones(n,1);   % multiplicador associado a restricao l-x
w=v;         % multiplicador associado a restricao x-u
parini=1;
r=parini;          % Parametro de penalidade inicial
%r=100;
%lamb = -776.3768;
%lamb =  -2.0628e+03;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Principal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timen=0;
tt=cputime; k=1;  aux=0;  % Time
%calculando P inversa b:
pinvb = P\b; pinva = P\a; rest0 = b'*xk0 - c; rest=rest0; xk=xk0; erro4=1;
while (k<100)&&(aux==0)
  
% %   bt=sqrt(r)*b;
% %   alfa=1/(1+bt'*bt);
% %   IP= eye(n)-(alfa*(bt*bt')); % inversa de P formula de Sherman-Morrisson
% %   xk=IP*(a-((lamb-r*c)*b)); % Solucao de grad(lag)=0

alfa=r/(1+r*b'*pinvb);
lar = r*c - lamb;

difproj = xk0 - (P*xk0 + (b'*xk0)*r*b - a - lar*b);
difproj=max(l,min(difproj,u));
erro3=norm(xk0-difproj);
erroold=erro3;
if flag>=10
   erroold=erro4; 
end

if (flag~=5)
    xk = pinva + (lar - alfa*(b'*pinva + lar*(b'*pinvb)))*pinvb;

    xk=max(l,min(xk,u));

    rest=b'*xk-c;
    %pause
end

 %comparar o lambcompara
    if flag==9
    matP=[]; sizm=1; reduxb=[]; reduxa=[];
    for i=1:n
       if (xk(i)<u(i))&&(xk(i)>l(i))
        matP(sizm,:) = P(i,:);
        reduxb(sizm,:) = b(i);
        reduxa(sizm,:) = a(i);
        sizm=sizm+1;
       end
    end
    if sizm>1 %ai a matP n ficou vazia n?...
       ovP = zeros(sizm,n+1);
       ovP(1:sizm-1,1:n) = matP;
       ovP(1:sizm-1,n+1) = reduxb;
       ovP(sizm,1:n) = b';
    end
    cvex = zeros(sizm,1);
    cvex(1:sizm-1) = reduxa;
    cvex(sizm) = c;
    inver = ovP'*ovP;
    solvvec = inver\(ovP'*cvex);
    lambcompara = solvvec(n+1);
    erro4 = norm(lamb-lambcompara);
    pareps = erro4/erroold;
    erroold = erro4;
    end
    %diflambs = lambcompara - lamb    
    %pause
    
    %lambcompdois, dessa vez com o mu
    if flag>=10
        veca = a - P*xk;
        %procurando as ativas
        ativ=0; bbar=-777*ones(2,1); vecabar = bbar; mu=zeros(n,1); indices = [];
        for i=1:n
           if (xk(i) == l(i))||(xk(i) == u(i))
              ativ = ativ+1;
              bbar(ativ) = b(i);
              vecabar(ativ) = veca(i);
              indices = [indices; i];
           end
        end
           %vendo qdo da o ativ e separando em casos
           if ativ == 0 %pto no interior
               lambcompdois = (b'*veca)/(b'*b);
               %erro4 = norm(-veca + lambcompdois*b);
               erro4 = norm([-veca + lambcompdois*b;b'*xk-c]);
               %pause
           end
           if ativ ==1
               bbar = [bbar(1)]; vecabar=[vecabar(1)];
           end
           if (ativ>=1)&&(ativ<n)
               %ai vc resolve o sistema
               %matriz = eye(ativ+1,ativ+1);
               %matriz(n,1:ativ) = bbar';
               %matriz(1:ativ,n) = bbar;
               %matriz(ativ+1,ativ+1) = b'*b;
               %vecabar = [vecabar; b'*veca];
               %meuvet = matriz\vecabar;
               
               %lambcompdois = meuvet(ativ+1);
               %mua = meuvet(1:ativ);
               %for i=1:n
               %   for j=1:length(indices)
               %       if i==indices(j)
               %           mu(i) = mua(j);
               %       end
               %   end
               %end
               
               %outro teste
               lambcompdois = (b'*veca - bbar'*vecabar(1:ativ))/((b'*b) - (bbar'*bbar));
               for i=1:n
                   if (xk(i)==l(i))
                       mu(i) = veca(i) - lambcompdois*b(i);
                       mu(i) = min(0,mu(i)); %para ficar -mu_l
                   end
                   if (xk(i)==u(i))
                       mu(i) = veca(i) - lambcompdois*b(i);
                       mu(i) = max(0,mu(i)); %para ficar mu_u
                   end
               end
               
               %erro4 = norm(-veca + mu + lambcompdois*b);
               erro4 = norm([-veca + mu + lambcompdois*b;b'*xk-c]);

           end
           if ativ==n
               %ai vc resolve o sistema
               matriz = eye(ativ+1,ativ+1);
               matriz(n,1:ativ) = bbar';
               matriz(1:ativ,n) = bbar;
               matriz(ativ+1,ativ+1) = b'*b;
               vecabar = [vecabar; b'*veca];
               meuvet = matriz\vecabar;
               
               lambcompdois = meuvet(ativ+1);
               mua = meuvet(1:ativ);
               for i=1:n
                  for j=1:length(indices)
                      if i==indices(j)
                          mu(i) = mua(j);
                      end
                  end
               end
               %erro4 = norm(-veca + mu + lambcompdois*b);
               erro4 = norm([-veca + mu + lambcompdois*b;b'*xk-c]);

           end %ativ n
    end %flag 10
    
    
    
    %if (abs(rest)>0.1*abs(rest0))
    %    r=2*r;
    %end   
    
    if (flag==2)&&(r>10000)
       aux=2; 
       %xk0 = quadprog(P,-a,[],[],b',c,l,u);
       %[xk0,fk0,~,out]=fmincon(@(x) x'*(P*x)-a'*x,xk0,[],[],b',[c],l,u);
    end
    if  (  (  (flag==3)&&(erro3>pareps)&&(flag<10)  )||(  (flag>=6)&&(erro3>pareps)&&(flag<10)  )||(flag==5)||(  (flag>=10)&&(erro4>pareps)  )  )&&(k>=2)
        %(  (flag==3)&&(r>100)  )||(flag==5)
        
       %nesse caso recomputo xk mas continuo o processo
       %xk = wright(P+r*(b*b'),-a+(lamb-r*c)*b,l,u,xk);
       
       if (flag==3)||(flag>=6)
        r=parini; 
       end

       options = optimoptions(@quadprog,'Algorithm','trust-region-reflective'); 
       xk = quadprog(P+r*(b*b'),-a+(lamb-r*c)*b ,[],[],[],[],l,u,xk,options);
       %xk = quadprog(P+r*(b*b'),-a+(lamb-r*c)*b ,[],[],[],[],l,u);
       xplot=xk(1); yplot=xk(2);
       %pause
       %xk = quadprog(P+r*(b*b'),-a+(lamb-r*c)*b ,[],[],[],[],l,u);
      
       if flag~=11
           flag=min(5,flag);
       else
           aux= 4;
       end
       rest=b'*xk-c;
       
       %xk0 = quadprog(P,-a,[],[],b',c,l,u);
       %[xk0,fk0,~,out]=fmincon(@(x) x'*(P*x)-a'*x,xk0,[],[],b',[c],l,u);
    end
    rest0 = rest;
        
    %r=n/(4*k+1);x0
   
    %r=abs(lamb-lamb0)/abs(rest0)
      k=k+1;
      
      %% Criterio de parada
      timen=timen+cputime-tt;
      
    erro1=norm(xk-xk0);
    erro2=abs(rest);
    difproj = xk - (P*xk + (b'*xk)*r*b - a - lar*b);
    difproj=max(l,min(difproj,u));
    erro3=norm(xk-difproj);
    %razao = norm(erro3/lamb);
    %norms = norm(xk0-xquadprog);
    
    
        %lamb0=lamb;
    lamb=lamb+r*rest;
    lar = r*c - lamb;
    
    
    if flag<6
      pareps = max(pardesc*pareps,prec);
    else
        if flag == 6
           pareps = max(erroold*0.9,prec);
        end
        if flag == 7
           pareps = max(erroold*0.5,prec);
        end
        if flag == 8
           pareps = max(erroold*0.1,prec);
        end
    end
     erroold=erro3;
    
    tt = cputime;
    if (erro1 <= prec && erro2 <= prec)&&(  ( erro3<=prec && flag==4  )||(flag~=4)  )&&(  ( erro4<=prec && flag>=9  )||(flag<9)  )
        aux=1;
    else
        if aux~=2
           xk0=xk; 
        else
            if aux~=4
              xk0 = quadprog(P,-a,[],[],b',c,l,u);
            else
                %xk0 = bi_section_fuel_v1(P,a,b,c,l,u);
                %xk0 = secant_fuel_v1(P,a,b,c,l,u);
                %xk0 = regula_falsi_fuel_v1(P,a,b,c,l,u);
                xk0 = regula_falsi_fuel_v1(P,a,b,c,l,u,lambcompdois);
            end
           %disp('linha 310');
           %[xk0,fk0,~,out]=fmincon(@(x) x'*(P*x)-a'*x,xk0,[],[],b',[c],l,u); 
        end
      
      %r=10;
    end
    %pause
end
epslon = erro3;
%difproj = xk - (P*xk + (b'*xk)*r*b - a - lar*b);
%difproj=max(l,min(difproj,u));
%erro3=norm(xk-difproj);
%epslon = erro3
%pause

%[xk0,fk0,~,out]=fmincon(@(x) x'*(P*x)-a'*x,xk0,[],[],b',[c],l,u);
%xk0 = quadprog(P,-a,[],[],b',c,l,u,xk0);

%plot(xplot,yplot,'.');


%lamb
%xk
%k


%xk=disp([xk(1); xk(2)])
%disp([erro1 erro2 timen])
%errot = norm(p.*xs - a + lambxs.ineqlin*b + lambxs.upper - lambxs.lower)
%disp([k]);

%comparando quadprog:
tic
xquadprog = quadprog(P,-a,[],[],b',c,l,u);
if length(xquadprog) == 0
   xquadprog = Inf*(ones(n,1)); 
end
timeprog=toc;
%[xquadprog,fquad]=fmincon(@(x) 0.5*x'*(P*x)-a'*x,ones(n,1),[],[],b',[c],l,u);
