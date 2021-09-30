purge

load planets.mat

strings = [64,59,56,52,47,44,40,35,32,28];

tsound = [];
notes = [];

for ii = 2:length(tout)-1
  for i = 1:Numparticles
%     test1 = sign(xpos_vec(i,ii)) ~= sign(xpos_vec(i,ii+1));
%     test2 = sign(ypos_vec(i,ii)) ~= sign(ypos_vec(i,ii+1));
%     if test1 || test2
%       tsound = [tsound;tout(ii)];
%       notes = [notes;strings(i)];
%     end
    for j = i:Numparticles
      if i ~= j
	distnow = sqrt((xpos_vec(i,ii)-xpos_vec(j,ii))^2+ ...
		       (ypos_vec(i,ii)-ypos_vec(j,ii))^2);
	distprev = sqrt((xpos_vec(i,ii-1)-xpos_vec(j,ii-1))^2+ ...
		       (ypos_vec(i,ii-1)-ypos_vec(j,ii-1))^2);
	distnext = sqrt((xpos_vec(i,ii+1)-xpos_vec(j,ii+1))^2+ ...
		       (ypos_vec(i,ii+1)-ypos_vec(j,ii+1))^2);
	
	if distnow < distprev && distnow < distnext
	  if distnow < 10
	    tsound = [tsound;tout(ii)];
	    notes = [notes;strings(j)];
	  end
	end
      end
    end
  end
end

plot(tsound,notes,'bo')
figure()
plot(xpos_vec',ypos_vec')
outvec = [tsound,notes];

save('notes.mat','outvec')