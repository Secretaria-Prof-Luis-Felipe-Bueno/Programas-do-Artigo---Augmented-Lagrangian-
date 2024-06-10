function x = xis_fuel_v1(L,P,a,b,l,u)
n = size(u);

x = -P\(L*b-a);
for j=1:n
    if x(j) < l(j)
        x(j) = l(j);
    else
        if x(j) > u(j)
            x(j) = u(j);
        end
    end
end
end