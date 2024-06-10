n = 500; m = 500;
for i=1:100
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matrizes1sobre',num2str(n),num2str(m),strmat);
    lp = rand(m,m);
    % grande=zeros(m,1);
    % for j=1:m grande(j)=j; end
    lp=(1/100)*(lp')*lp + eye(n,n);
    P = zeros(n,n);
    %for i=1:m P(i,i) = rand(1); end
    P(1:m,1:m) = lp;
    a = rand(n,1);
    save(mystr,'P','a');    
end
