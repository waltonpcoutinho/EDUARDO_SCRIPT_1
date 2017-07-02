function [fbest pbest fconv] = GA_clu(data_path, qtd, k, dim, TimeToRun)

   format long g

   tic
   t = toc;

   range=[0 0 qtd-1 dim-1];
   data = dlmread(data_path,',',range);

   % Captura dos valores minimos e maximos para cada dimensão de acordo com os limites dos dados de entrada
   for i = 1:dim;
       [minimo,imin]=min(data(:,i));
       [maximo,imax]=max(data(:,i));
       l(i)=minimo;
       u(i)=maximo;
   end;

   % define as mesmas fronteiras para cada centro
   for i = 1:k-1
      l(i*dim+1:i*dim+dim) = l(1:dim);
      u(i*dim+1:i*dim+dim) = u(1:dim);
   end

   %nGen = 100;
   N = 100; % population size
   n = k*dim; % tamanho de cada cromossomo
   P = zeros(N, n); % cria populacao
   pcross = 0.8;
   pmut = 0.001;

   [dummy, ia, ic] = unique(data,'rows');
   nCand = size(ia);
   nCand = nCand(1);

   nimprov = 1;

   fbest = inf;
   for i=1:N % inicializa populacao aleatoriamente

      rObjs = randperm(nCand,k);
      rObjs = ia(rObjs);

   	  for c=0:k-1
   	  	P(i, (c*dim+1):(c*dim+dim)) = data(rObjs(c+1),:);
   	  end

   	  [fitness(i) P(i,:)] = fitnessCalc(k, qtd, dim, data, P(i,:));
   	  if fitness(i) < fbest
   	  	 fbest = fitness(i);
         pbest = P(i,:);
         fconv(nimprov,1) = toc;
         fconv(nimprov,2) = fbest;
         nimprov = nimprov + 1;
   	  end

      if ttest(TimeToRun); return; end; % fim da execução?

   end

   %for i=1:nGen
   while t < TimeToRun

      selected = zeros(N);
      rouletteCirc = sum(fitness);
      fitnessAux = fitness;

      for j=1:(pcross*N/2); % selection and crossover
         
         p1 = roulette_wheel_selection(rouletteCirc, fitnessAux, N);
         while selected(p1) == 1
            p1 = roulette_wheel_selection(rouletteCirc, fitnessAux, N);
         end
         selected(p1) = 1;
         p2 = roulette_wheel_selection(rouletteCirc, fitnessAux, N);
         while selected(p2) == 1
            p2 = roulette_wheel_selection(rouletteCirc, fitnessAux, N);
         end
         selected(p2) = 1;
        
         [P(p1,:) P(p2,:)] = crossover(P(p1,:), P(p2,:),n); 
         P(p1,:) = viability(P(p1,:), l, u, k, dim);
         P(p2,:) = viability(P(p2,:), l, u, k, dim);

         [fitness(p1) P(p1,:)] = fitnessCalc(k, qtd, dim, data, P(p1,:));
         [fitness(p2) P(p2,:)] = fitnessCalc(k, qtd, dim, data, P(p2,:));

         if(fitness(p1) < fbest) 
           fbest = fitness(p1);
           pbest = P(p1,:);
           fconv(nimprov,1) = toc;
           fconv(nimprov,2) = fbest;
           nimprov = nimprov + 1;
         end

         if(fitness(p2) < fbest) 
           fbest = fitness(p2);
           pbest = P(p2,:);
           fconv(nimprov,1) = toc;
           fconv(nimprov,2) = fbest;
           nimprov = nimprov + 1;
         end

         if ttest(TimeToRun); return; end; % fim da execução?

      end

      for j=1:N % mutation

         aux = P(j,:);

         for g=1:n
           
           if (rand > pmut) continue; end;
           
           delta = rand;
           
           if P(j,g) ~= 0

              sign = rand;
              if (sign > 0.5)
                P(j,g) = P(j,g) + 2 * delta * P(j,g);
              else
                P(j,g) = P(j,g) - 2 * delta * P(j,g);
              end

           else 

              sign = rand;
              if (sign > 0.5)
                P(j,g) = P(j,g) + 2 * delta;
              else
                P(j,g) = P(j,g) - 2 * delta;
              end

           end
            
         end

         [fitness(j) P(j,:)] = fitnessCalc(k, qtd, dim, data, P(j,:));

         if(fitness(j) < fbest) 
           fbest = fitness(j);
           pbest = P(j,:);
           fconv(nimprov,1) = toc;
           fconv(nimprov,2) = fbest;
           nimprov = nimprov + 1;
         end

         if ttest(TimeToRun); return; end; % fim da execução?

      end

      t = toc;

   end % main loop


end

function C = viability(P, l, u, k, d)

  dim = 1;
  for i = 1:k;
     for j = 1:d;
        if P(dim) < l(j)
          P(dim) = l(j);

        elseif P(dim) > u(j)
          P(dim) = u(j);
        end
        dim = dim +1;
     end
  end

  C = P;

end

function [C1 C2] = crossover(P1, P2, n) 

  cp = randi([2 n]); % 1-cut point
  C1 = P1;
  C2 = P2;

  for i=cp:n
     C2(i) = P1(i);
     C1(i) = P2(i);
  end

end

function p = roulette_wheel_selection(rouletteCirc, fitness, N)

  r = randi([0 uint64(rouletteCirc)]);

  i=1;
  sum = fitness(i);
  while sum < r
    i = i+1;
    if(i > N) 
      i = i-1;
      break; 
    end;
    sum = sum + fitness(i);
  end

  p = i;

end

function [value centroids] = fitnessCalc(k, qtd, d, data, p)

   soma_dist_c = zeros(k,1); cont_c = zeros(k,1);
   soma_total = 0;  soma_d = zeros(k,d);
   cont_pontos_cent = zeros(1,k);
   sum_val_cent = zeros(k,d);

   for obj=1:qtd;
      dim = 0;
      for c = 1:k;
          soma_c = 0;
          for j = 1:d;
              dim = dim + 1;
              soma_c = soma_c + (data(obj,j) - p(dim))^2;
          end;
          dist_c(c) = sqrt(soma_c);
      end;
      [min_dist,imin_dist] = min(dist_c);
      cent_ponto(obj) = imin_dist;
      cont_pontos_cent(1,imin_dist) = cont_pontos_cent(1,imin_dist)+1;

      for j=1:d
      	sum_val_cent(imin_dist,j) = sum_val_cent(imin_dist,j) + data(obj,j);
      end

      soma_total = soma_total + min_dist;
   end;

   if any(cont_pontos_cent == 0)
      value = soma_total;
      centroids = p;
      return;
   end

   centroids = zeros(1,k*d);
   dim = 0;
   for c=1:k
     for j=1:d
       dim = dim+1;
       centroids(dim) = sum_val_cent(c,j)/cont_pontos_cent(c);
     end
   end

   value = soma_total;

end

function flag = ttest(timeout)
  ctime = toc;
  if ctime > timeout
    flag = true;
  else
    flag = false;
  end
end