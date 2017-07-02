function [Jb Ab fconv] = TS_clu(data_path, m, c, d, TimeToRun)

   format long g % define formato de saida

   tic % salva tempo atual
   runtime = toc; % tempo de execucao desde o ultimo tic

   range=[0 0 m-1 d-1]; % define limites dos dados [x1 y1 x2 y2] para leitura do arquivo
   x = dlmread(data_path, ',', range); % leitura da base de dados
   P = 0.99; % threshold probability
   MTLS = m; % tabu list max size
   NTS = 20; % number of neighbors generated
   ITMAX = 2000;
   TLL = 0; % tabu list current size
   TL = zeros(1,1); % tabu list 

   Ac = randi([1 c],[1 m]); % current solution
   Ab = Ac; % best solution is current solution
   Jc = fo(c,m,d,x,Ac); % objective function of current solution
   Jb = Jc; % best f.o is current f.o
   At = zeros(NTS, m); % trial solutions
   Jt = zeros(1,NTS); % objective function of the trial solutions

   nimprov = 1;
   %time_over = false;
   %for k=1:ITMAX
   while runtime < TimeToRun

      for t=1:NTS % generates NTS neighbors of Ac
        At(t,:) = generateNeighbor(Ac, m, c, P);
        Jt(t) = fo(c, m, d, x, At(t,:));
        runtime = toc;
        % if runtime >= TimeToRun
        % 	time_over = true;
        % 	break
        % end
        if ttest(TimeToRun); return; end; % fim da execução?

      end

      % if time_over
      % 	break
      % end

      [Jt I] = sort(Jt);

      all_TL = 1;
      for t=1:NTS
        if not(any(ismember(TL, Jt(t))))
          all_TL = 0;
          Ac = At(I(t),:);
          Jc = Jt(t);        
          break              
        end

        if ttest(TimeToRun); return; end; % fim da execução?

      end

      if all_TL == 0
         TLL = TLL + 1;
         if TLL == MTLS + 1
            TL(TLL) = Jc;
            TL(1) = [];
            TLL = TLL - 1;
         else
            TL(TLL) = Jc; % add solution in the tabu list
         end
         
         if Jc < Jb % updates best clustering
            Jb = Jc;
            Ab = Ac;
            fconv(nimprov,1) = toc;
            fconv(nimprov,2) = Jb;
            nimprov = nimprov + 1;
         end
      end

      runtime = toc;

   end

end

function At = generateNeighbor(Ac, m, c, P) 

  for o=1:m
    R = rand;
    if R < P
       At(o) = Ac(o);
    else
       cand = [1:c]; 
       cand(cand==Ac(o)) = [];
       id = randi([1 c-1]);
       At(o) = cand(id);
    end
  end

end

% calculates objective function
function [value cent] = fo(c, m, d, x, A)

  qtd_cent = zeros(1,c);
  sum_clust = zeros(c,d); % soma dos valores dos pontos em cada dimensao, separados por cluster 

  % calcula qtd_cent e sum_clust
  for o=1:m
     qtd_cent(A(o)) = qtd_cent(A(o)) + 1;
     for j=1:d 
        sum_clust(A(o),j) = sum_clust(A(o),j) + x(o,j);
     end
  end

  % centroids = zeros(c,d);
  % for i=1:c % calcula centroides
  %    if (qtd_cent(i) == 0)
  %      value = inf;
  %      return; 
  %    end % solucao inviavel pois tem cluster vazio
  %    centroids(i,:) = sum_clust(i,:)/qtd_cent(i);
  % end
  
  % iClusterDist = 0;  % sum of intracluster distances
  % % calcula soma das distancias intracluster
  % for o=1:m
  %   cent_dist = 0;
  %   for j=1:d
  %       cent_dist = cent_dist  + (x(o,j) - centroids(A(o),j))^2;
  %   end
  %   cent_dist = sqrt(cent_dist); % euclidian distance to the centroid
  %   iClusterDist = iClusterDist + cent_dist;
  % end

  % value = iClusterDist;

  cent = zeros(1,c*d);

  dim = 1;
  for i=1:c

	if (qtd_cent(i) == 0)
		value = inf;
		return 
	end 

    for j=1:d
      cent(dim) = sum_clust(i,j)/qtd_cent(i);
      dim = dim + 1;
    end

  end

  value = calcObjFunc(c, m, d, x, cent);

end

function value = calcObjFunc(k, l, d, data, x) 

   soma_dist_c = zeros(k,1); cont_c = zeros(k,1);
   soma_total = 0;  soma_d = zeros(k,d); 
   cont_pontos_cent = zeros(k);
   for obj=1:l;
      dim = 0;
      for c = 1:k;
          soma_c = 0;
          for j = 1:d;
              dim = dim + 1;
              soma_c = soma_c + (data(obj,j) - x(dim))^2;
          end;    
          dist_c(c) = sqrt(soma_c);
      end; 
      [min_dist,imin_dist] = min(dist_c);
      cent_ponto(obj) = imin_dist;
      cont_pontos_cent(imin_dist) = cont_pontos_cent(imin_dist)+1;
      soma_total = soma_total + min_dist;
   end;
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