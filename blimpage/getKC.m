purge
global K C

costmin = 1e20;
step = 0.05;
for K = 5:step:7
    K
    for C = 5:step:7
        cost = Contact(0);
        if cost < costmin
            costmin = cost;
            Kmin = K;
            Cmin = C;
        end
    end
end
disp('Solution')
K = Kmin
C = Cmin
disp('Diameter in Inches')
b = 38
Contact(1);
