function [fgbest gbest fconv] = PSO_clu(data_path, qtd, k, dim, TimeToRun)
   
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

   n = k*dim; % dimensao do espaço de busca
   N=40; % tamanho do enxame de partículas
   %n_iter = n*15; % número de iterações do PSO
   w_max = 0.4; % peso inercial máximo
   w_min = 0.1; % peso inercial mínimo
   c1 = 1.5; c2 = 1.5; % constantes de aceleração do componente cognitivo e social, respectivamente
   %w = linspace(w_max, w_min, n_iter); % gera vetor com n_iter pesos inerciais igualmente espaçados decrescendo de w_max para w_min  
   p = zeros(N, n); % cria enxame de partículas p_i com os valores de centros candidatos
   v = zeros(N, n); % vetor velocidade das partículas inicializado com zero

   for i=1:N
      for d=1:n
         p(i,d) = l(d) + rand * (u(d) - l(d)); % inicializa enxame aleatoriamente com limiares em tons de cinza entre o intervalo presente na imagem
         v(i,d) = l(d) + rand * (u(d) - l(d)); % inicializa velocidades aleatoriamente
      end
   end
   
   for i=1:N
      fpbest(i) = calcObjFunc(k, i, qtd, dim, data, p); % calcula função objetivo       
      pbest(i,:) = p(i,:);
   end

   [fgbest igbest] = min(fpbest); % captura o valor da função objetivo e o indice da particula gbest 
   gbest = p(igbest,:); % salva particula gbest
  
   ctime = toc;
   firstImp = true;

   nimprov = 1;
   % laço principal de iterações do PSO
   % cada iteração consiste em uma movimentação do enxame
   it = 1;
   while t < TimeToRun
   % for t=1:n_iter

      w = w_max - (t/TimeToRun * (w_max - w_min));

      for i=1:N

         for j=1:n
            %v(i,j) = w(t) * v(i,j) + c1*rand*(pbest(i,j) - p(i,j)) + c2*rand*(gbest(j) - p(i,j)); % calcula componente j do vetor velocidade
            v(i,j) = w * v(i,j) + c1*rand*(pbest(i,j) - p(i,j)) + c2*rand*(gbest(j) - p(i,j)); % calcula componente j do vetor velocidade
            p(i,j) = p(i,j) + v(i,j); % movimenta particula na componente j dado o j do vetor velocidade calculado
            % verifica se a particula desreipeitou os limites de tons de cinza e faz o recuo se necessário
            
            if p(i,j) < l(j) 
               p(i,j) = l(j);
            elseif p(i,j) > u(j)
               p(i,j) = u(j);
            end
         end
         
         fo = calcObjFunc(k, i, qtd, dim, data, p); % calcula f.o da nova posicao da particula
         
         if fo < fpbest(i) % atualiza o pbest se necessario
            fpbest(i) = fo;
            pbest(i,:) = p(i,:);
         end
         if fo < fgbest % atualiza o gbest se necessário
            gbest = p(i,:);
            igbest = i;
            fgbest = fo;
            if toc - ctime > 5 || firstImp
              fconv(nimprov,1) = toc;
              fconv(nimprov,2) = fgbest;
              nimprov = nimprov + 1;
              ctime = toc;
              firstImp = false;
            end
         end   

         if ttest(TimeToRun); return; end; % fim da execução?

      end

      t = toc;
      it = it +1;

   end % main loop
   
   %T = sort(gbest);
end

function value = calcObjFunc(k, part, l, d, data, x) 

   soma_dist_c = zeros(k,1); cont_c = zeros(k,1);
   soma_total = 0;  soma_d = zeros(k,d); 
   cont_pontos_cent = zeros(k);
   for obj=1:l;
      dim = 0;
      for c = 1:k;
          soma_c = 0;
          for j = 1:d;
              dim = dim + 1;
              soma_c = soma_c + (data(obj,j) - x(part,dim))^2;
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
