r1 = ceil((N).*rand(1,1));
while r1 >= N
  r1 = ceil((N).*rand(1,1));
end
r2 = ceil((N).*rand(1,1));
while r2 <= r1
  r2 = ceil((N).*rand(1,1));
end
cross = [xchild(2,r1:r2);xchild(1,r1:r2)];
xchild(:,r1:r2) = cross;