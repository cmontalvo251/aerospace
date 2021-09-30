%%Genetic Alogrith for a glider
clear
clc
close all

%%Generate Initial Parameters
parameters

%%Create Initial Random Population
xbin = randi(10,[popsize,N])

disp('Initial Population')
disp(xbin)
disp(' ')

iteration = -1;

while iteration < maxiteration

   Fi = zeros(popsize,1);  
   for ii = 1:popsize %%Calculate Fitnesses
      Fi(ii) = fitness(xbin(ii,:));
   end
   
   %%Search for Elitist1
   maxloc = max(find(Fi == max(Fi)));
   elite1 = xbin(maxloc,:);
   Fi2 = Fi;Fi2(maxloc)=0;
   %%Search for Elitist2
   maxloc2 = max(find(Fi2 == max(Fi2)));
   elite2 = xbin(maxloc2,:);
   
   disp('Current Elitist')
   disp(elite1);
   disp(elite2);
   disp(' ')
   
   %%Construct Roulette Wheel
   Fsum = sum(Fi);
   Pwheel =(Fi./Fsum);
   Pmin = min(Pwheel);
   Pwheel = cumsum(Pwheel);
   
   %%%Create Next Generation
   xgen = xbin;parents = zeros(1,2);
   for i = 1:popsize/2-1
     
     %%pick 2 different parents
      while parents(1) == parents(2) 
         s1 = rand;s2 = rand;
         while s1 > max(Pwheel) || s2 > max(Pwheel)
            s1 = rand;s2 = rand; %%dont overshoot Pwheel
         end
         c = 0;parents = zeros(1,2);
         while parents(1) == 0 || parents(2) == 0
            c = c + 1;
            if  parents(1) == 0 && s1 < Pwheel(c)
               parents(1) = c;
            end
            if  parents(2) == 0 && s2 < Pwheel(c)
               parents(2) = c;
            end
         end
      end
      Pwheel(parents) = [];
      xchild = xbin(parents,:);
      xbin(parents,:) = [];
      parents = zeros(1,2);
      
      %%%Check For Possibility of Crossover between parents
      if rand < Pcross 
	crossover
      end
      
      %%Check for Possibility of Mutation
      for k = 1:2
         if rand < Pmutation 
	   mutate
         end
      end
      xgen([2*i-1 2*i],:) = xchild;
   end
   xgen([end-1,end],:) = [elite1;elite2];
   xbin = xgen;
   
   disp('Next Generation')
   disp(xgen)
   disp(' ')

   pause
end
