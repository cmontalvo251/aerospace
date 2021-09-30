r = ceil((N).*rand(1,1));	    
dev = (delMutate/100)*xchild(k,r);
xchild(k,r) = xchild(k,r) + dev*randn;