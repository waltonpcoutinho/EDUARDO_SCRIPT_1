% Continuous Grasp (C-GRASP) from "GLOBAL OPTIMIZATION BY CONTINUOUS GRASP", M. J. HIRSCH, C. N. MENESES, P. M. PARDALOS, AND M.G.C. RESENDE, 2006
% n - number of dimensions
% l - lower bound vector
% u - upper bound vector, then l <= x <= u

function [f_best x_best fconv] = CGRASP_clu(data_path, qtd, k, dim, TimeToRun, MaxDirToTry, delta1, delta2)
   
   format long g
   
   tic;
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
   
   h_max = max(means);
   h_min = min(means)/10000;
   I_s = 100; % número de tentativas da inicialização

   % define as mesmas fronteiras para cada centroide 
   for i = 1:k-1
      l(i*dim+1:i*dim+dim) = l(1:dim);
      u(i*dim+1:i*dim+dim) = u(1:dim);  
   end

   n = k*dim; % dimensao do espaço de busca

   f_best = inf; 

   [dummy, ia, ic] = unique(data,'rows');
   nCand = size(ia);
   nCand = nCand(1);

   nimprov = 1; % auxiliar para grafico de convergencia
   
   while t < TimeToRun

      improve = false; % para verificar se f_best melhorou na construcao
      firstimprove = false;

      % gera aleatoriamente ponto inicial entre l e u que definem as fronteiras da região
      % for d = 1 : n 
      %    x(d) = l(d) + rand * (u(d) - l(d));
      % end

      fbest_x = inf;
      for it=1:I_s % inicializa populacao aleatoriamente

         rObjs = randperm(nCand,k);
         rObjs = ia(rObjs);

         for c=0:k-1
         	x(1, (c*dim+1):(c*dim+dim)) = data(rObjs(c+1),:);
         end

		     [fo x] = fo_cent(k, qtd, dim, data, x);

         %fo = f(x, n, k, data, qtd);
         if fo < fbest_x
           fbest_x = fo;
           xbest = x;
           if fo < f_best
               f_best = fo;
               x_best = x;
               fconv(nimprov,1) = toc;
               fconv(nimprov,2) = f_best;
               nimprov = nimprov + 1;
           end
         end

         if ttest(TimeToRun); return; end; % fim da execução

      end

      x = xbest;
      if f_best > fbest_x
         f_best = fbest_x;
         x_best = x;
         fconv(nimprov,1) = toc;
         fconv(nimprov,2) = f_best;
         nimprov = nimprov + 1;
      end

      %[fo x] = fo_cent(k, qtd, dim, data, x);

      h = h_max; % controla a discretização do espaço de busca
      
      f_ant = inf; % retem a solução obtida na iteração anterior caso nao seja a primeira iteração
      jump = false; % define se a iteração nao é mais promissora e pula para a próxima
      
      % Multi-start sobre a construção para h de hmax a hmin
      while h > h_min
         
         [x ImprC] = ConstructGreedyRandomizedImproved(x, n, h, l, u, k, data, qtd, TimeToRun);

         fo = f(x, n, k, data, qtd);

         if fo  < f_best
            x_best = x; f_best = fo;
            improve = true;
            fconv(nimprov,1) = toc;
            fconv(nimprov,2) = f_best;
            nimprov = nimprov + 1;
         end

         if f_ant == inf
            f_ant = fo;
         elseif ((f_ant-fo == 0) & (firstimprove == true)) | (f_ant-fo ~= 0) % so usa filtro se houve a primeira melhora 
            firstimprove = true;
            gap = (f_ant-fo)*100/f_ant;
            if gap < delta1
               gap_b = (fo-f_best)*100/fo;
               if gap_b > delta2
                  jump = true;
                  break;
               end;
            end;
            f_ant = fo;
         end;
         
         h = h/2; % aumenta densidade da grade explorada

         if ttest(TimeToRun); return; end; % fim da execução?

      end
   
      if t > TimeToRun
         break;
      end;
      
      if jump == true | improve == false
         continue;
      end;
      
      H = 1;
      while H > (h_min/10)

         [x fconv nimprov] = LocalSearch(x_best, n, H , l, u, MaxDirToTry, k, data, qtd, fconv, nimprov, TimeToRun);
         fo = f(x, n, k, data, qtd);

         if fo < f_best
            x_best = x; f_best = fo;
            fconv(nimprov, 1) = toc;
            fconv(nimprov, 2) = f_best;
            nimprov = nimprov + 1;
         end
         
         H = H/2;

         if ttest(TimeToRun); return; end; % fim da execução?

      end

      t = toc; % current elapsed time

   end

end


function [value centroids] = fo_cent(k, qtd, d, data, p)

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