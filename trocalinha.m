function [y] = trocalinha(stru,i,j)

teststruc = stru;
[n,~] = size(stru);
n=n-1;
for k=1:n
    teststruc{k,i}=stru{k,j};
end
for k=1:n
    teststruc{k,j}=stru{k,i};
end
y = teststruc;

end

