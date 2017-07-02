% Continuous Grasp (C-GRASP) from "GLOBAL OPTIMIZATION BY CONTINUOUS GRASP", M. J. HIRSCH, C. N. MENESES, P. M. PARDALOS, AND M.G.C. RESENDE, 2006
% n - number of dimensions
% l - lower bound vector
% u - upper bound vector, then l <= x <= u

function [f_best x_best fconv] = CGRASP_2007(data_path, qtd, k, dim, TimeToRun, MaxIters, MaxNumIterNoImprov, MaxDirToTry, alpha)
   
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
       means(i) = mean(data(:,i));
   end;
   
   % define as mesmas fronteiras para cada centroide 
   for i = 1:k-1
      l(i*dim+1:i*dim+dim) = l(1:dim);
      u(i*dim+1:i*dim+dim) = u(1:dim);  
   end

   n = k*dim; % dimensao do espaço de busca

   f_best = inf; 
   
   h_max = 1;
   nimprov = 1;
   fconv(nimprov,1) = 0;

   while t < TimeToRun
      
      % gera aleatoriamente ponto inicial entre l e u que definem as fronteiras da região
      for d = 1 : n 
         x(d) = l(d) + rand * (u(d) - l(d));
      end
      h = h_max; NumIterNoImprov = 0;% controla a discretização do espaço de busca
      
      for Iter = 1 : MaxIters
         x = ConstructGreedyRandomized(x, n, h, l, u, alpha, k, data, qtd, TimeToRun);
         [x fconv nimprov] = LocalSearch(x, n, h , l, u, MaxDirToTry, k, data, qtd,fconv, nimprov, TimeToRun);
         
         fo = f(x, n, k, data, qtd);
         
         if fo  < f_best
            x_best = x; f_best = fo;
            fconv(nimprov,1) = toc;
            fconv(nimprov,2) = f_best;
            nimprov = nimprov + 1;
         else
            NumIterNoImprov = NumIterNoImprov + 1;   
         end
         
         if NumIterNoImprov >= MaxNumIterNoImprov
            h = h/2; NumIterNoImprov = 0; % torna grade o dobro mais densa         
         end
         
         if ttest(TimeToRun); return; end; % fim da execução?

      end;
      
   end
   
end

function flag = ttest(timeout)
  ctime = toc;
  if ctime > timeout
    flag = true;
  else
    flag = false;
  end
end
