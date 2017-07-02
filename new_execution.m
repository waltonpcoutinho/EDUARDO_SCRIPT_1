%function new_execution(method, dataset_id)
function new_execution(method, dataset_id)
    
   format long g;
   
   dpaths = {'datasets/ionosphere.data','datasets/wdbc.data', 'datasets/banknote.data', 'datasets/yeast.data', ...
             'datasets/dim032.data', 'datasets/dim064.data', 'datasets/a1.data', 'datasets/D31.data', ...
             'datasets/s1.data', 'datasets/unbalance.data', 'datasets/abalone.data', 'datasets/banknote.data', ... 
             'datasets/artset3.data', 'datasets/artset4.data', 'datasets/artset5.data', 'datasets/artset6.data',...
             'datasets/cmc.data', 'datasets/wine.data', 'datasets/vowel.data', 'datasets/breast-cancer-wisconsin.data',...,
              'datasets/iris.data','datasets/glass.data', 'datasets/crudeoil.data', 'datasets/thyroid.data',...,
              'datasets/artset1.data', 'datasets/artset2.data' }; 
   dnames = {'ionosphere', 'wdbc', 'banknote', 'yeast', 'dim32', 'dim64', 'a1', 'D31', 's1', 'unbalance', 'abalone', ...
             'banknote2','artset3', 'artset4', 'artset5', 'artset6','cmc', 'wine', 'vowel', 'cancer', 'iris', 'glass', ...,
             'crudeoil', 'thyroid', 'artset1', 'artset2'};
   
   lengths = {351, 569, 1372, 1484, 1024, 1024, 3000, 3100, 5000, 6500, 4177, 1372, 1500, 1750, 2000, 2500,...,
   				1473, 178, 871, 683, 150, 214, 56, 215, 250, 600}; % datasets length
   k = {2, 2, 10, 10, 16, 16, 20, 31, 15, 8, 8, 4, 3, 4, 5, 6, 3, 3, 6, 2, 3, 6, 3, 2, 5, 4}; % k clusters task
   dim = {34, 30, 4, 8, 32, 64, 2, 2, 2, 2, 3, 4, 2, 3, 4, 5, 9, 13, 3 ,9 ,4, 9, 5, 5, 3, 2}; % data dimension

   iter_num = 50;
   timeout = {1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,200,300,300,200, 30, 600, 10, 60,...,
   			  5, 5};

   I_s = 100; maxdirtotry = 30; delta1 = 0.1; delta2 = 0.5; %%% CGRASP-Clu parameters %%%
   alfa = 0.4; maxnumiternoimprov = 20; maxiters = 200; %%% CGRASP-2007 parameters %%%

   id = dataset_id;

   fo = zeros(1, iter_num); % fbest para 50 execucoes
   icc = zeros(1, iter_num); % icc para 50 execuções do CGRASP-clu
   bestcent = zeros( 1, k{id} * dim{id} ); % mantem melhor centroid encontrado nas exec. do CGRASP-clu
   ftimes = cell(1, iter_num); % curvas de convergencia para as 50 execuções
   fbest = inf; % melhor valor para a f.o. encontrado
   
   for it = 1 : iter_num
     	
     fval = 0; x = []; fconv = [];
      
     switch method
        case 'CGRASP_clu'
           [fval x fconv] = CGRASP_clu(dpaths{id}, lengths{id}, k{id}, dim{id}, timeout{id}, maxdirtotry, delta1, delta2);
        case 'CGRASP_2007'
           [fval x fconv] = CGRASP_2007(dpaths{id}, lengths{id}, k{id}, dim{id}, timeout{id}, maxiters, maxnumiternoimprov, maxdirtotry, alfa); 
        case 'GA_clu'
           [fval x fconv] = GA_clu(dpaths{id}, lengths{id}, k{id}, dim{id}, timeout{id}); 
        case 'TS_clu'
           [fval IDX fconv] = TS_clu(dpaths{id}, lengths{id}, k{id}, dim{id}, timeout{id});
           IDX = IDX';
        case 'PSO_clu'
           [fval x fconv] = PSO_clu(dpaths{id}, lengths{id}, k{id}, dim{id}, timeout{id});
        case 'kmeans'
           [fval IDX fconv] = kmeans(dpaths{id}, lengths{id}, dim{id}, k{id});
     end 
     
     fo(1, it) = fval; ftimes{1, it} = fconv; % registra f.o obtida na execução e curva de conv.
     
     if strcmp(method, 'TS_clu') || strcmp(method, 'kmeans')
     	icc(1, it) = ICC_idx(dpaths{id}, lengths{id}, dim{id}, IDX); % calcula ICC do resultado
     else
     	icc(1, it) = ICC(dpaths{id}, lengths{id}, k{id}, dim{id}, x); % calcula ICC do resultado
 	 end

     if strcmp(method, 'CGRASP_clu')
       if fval < fbest; bestcent(:) = x(:); fbest = fval; end;  % atualiza bestcent
     end 

   end;

   avgFconv = getAvgFconv( ftimes(1,:), timeout ); % gera grafico de convergencia medio
   string = strcat( 'results_v2/', method, '/', method ,'_', dnames{id} ); string = strcat(string, '.dat');
   dlmwrite(string, avgFconv, 'delimiter', '\t', 'precision', '%.3f'); % salva em arquivo a curva de convergencia media
	
   fileName = strcat( 'results_v2/', method, '/results_', dnames{id} );
   if strcmp(method, 'CGRASP_clu')
	   save(fileName, 'fo', 'bestcent', 'icc'); % salva em arquivo os resultados numéricos
   else
   	   save(fileName, 'fo', 'icc');
   end  
   
   exit();

end 

% calculates average convergence curve  
function avg_fconv = getAvgFconv(fconv, runtime) 
    
    [n niter] = size(fconv);
    avg_fconv = cell(1,6);
    
    for i=1:n
        
        % lista tempos envolvidos em todas as execucoes
        times = 0; count = 1;
        for t=1:niter
           [l c] = size(fconv{i,t});
           for k=1:l            
              times(1,count) = fconv{i,t}(k,1);  count = count+1;
           end
        end
        
        times = sort(times); % ordena tempos
        avg_fconv_i = zeros(1,1);
         
        % para cada tempo calcula a media entre todas as iteracoes 
        for c=1:count-1
           time = times(1,c);
           fos(:,:) = 0;
           for t=1:niter
              k = 1; fo = inf;
              [lin col] = size(fconv{i,t});
              while (k <= lin) & (fconv{i,t}(k,1) <= time)  
                fo = fconv{i,t}(k,2); k = k+1; 
              end
              fos(1,t) = fo;
           end
           avg_fconv_i(c,1) = time;
           avg_fconv_i(c,2) = mean(fos);
        end

        [lin col] = size(avg_fconv_i);
        avg_fconv_i(lin+1,1) = runtime{i};
        avg_fconv_i(lin+1,2) = avg_fconv_i(lin,2);
        
        aux = [0 avg_fconv_i(1,2)];
        for x=1:lin+1
           for y=1:col
              aux(x+1,y) = avg_fconv_i(x,y);  
           end
        end
        avg_fconv_i = aux;
        avg_fconv{1,i} = avg_fconv_i;
    end    

end

% calculates the Intracluster Correlation Coefficient by receive the data and the cluster centers
function icc = ICC(data_path, N, k, d, x)
  
  range=[0 0 N-1 d-1];
  data = dlmread(data_path,',',range);

  cdata = cell(1, k); % cada celula possui matriz com os pontos do cluster 

  for o = 1:N % para cada objeto, add linha em cdata para a matriz do seu cluster
    
    dim = 0;
    mi = inf; mi_id = -1;
    for c = 1:k;
      soma_c = 0;
      for j = 1:d
          dim = dim + 1;
          soma_c = soma_c + (data(o,j) - x(dim))^2;
      end    
      dist_c = soma_c;
      if dist_c < mi
        mi = dist_c;
        mi_id = c;
      end
    end 

    nextID = size( cdata{mi_id}, 1 ) + 1; % id da ultima linha na matriz do cluster mais um (nova linha)
    cdata{mi_id}(nextID, :) = data(o,:); % add ponto (linha) a matriz do seu cluster

  end

 
  w = zeros(1, k);  % probabilidade de cada cluster 
  m = zeros(k, d); % media de cada cluster (centroide)
  for i = 1:k
    w(1,i) = size(cdata{i}, 1)/N;
    m(i,:) = mean(cdata{i});
  end

  v_w = 0; % variancia intraclasse (within class variance)
  for i = 1:k
    v_i = 0;
    for j = 1:size( cdata{i}, 1 )
      ssq = 0; % sum of square error
      for l = 1:d
        ssq = ssq + ( cdata{i}(j,l) - m(i,l) )^2; 
      end
      v_i = v_i + ssq;
    end
    v_w = v_w + w(1,i) * v_i;
  end

  m_t = 0; % media total
  for i = 1:N
    m_t = mean(data);
  end

  v_t = 0; % total variance
  
  for i = 1:N % total variance
    ssq = 0; % sum of square error
    for j = 1:d
      ssq = ssq + ( data(i,j) - m_t(j) )^2; 
    end
    v_t = v_t + ssq;
  end

  %{
  v_b = 0;
  
  for i = 1:k
    ssq = 0;
    for j = 1:d
      ssq = ssq + ( m(i,j) - m_t(j) )^2;
    end
    v_b = v_b + w(i) * ssq; 
  end
  %}

  v_b = v_t - v_w;
  icc = v_b/(v_b + v_w);

end

% calculates the Intracluster Correlation Coefficient by receive the data and the cluster centers
function iccv = ICC_idx(data_path, N, d, idx)
  
  range=[0 0 N-1 d-1];
  data = dlmread(data_path,',',range);

  idx = idx+1;
  u = sort(unique(idx));
  k = size(u,1);
  cdata = cell(1, k); % cada celula possui matriz com os pontos do cluster 

  for o = 1:N % para cada objeto, add linha em cdata para a matriz do seu cluster
    cId = find(u == idx(o,1));
    nextID = size( cdata{cId(1,1)}, 1 ) + 1; % id da ultima linha na matriz do cluster mais um (nova linha)
    cdata{cId(1,1)}(nextID, :) = data(o,:); % add ponto (linha) a matriz do seu cluster
  end
 
  w = zeros(1, k);  % probabilidade de cada cluster 
  m = zeros(k, d); % media de cada cluster (centroide)
  for i = 1:k
    w(1,i) = size(cdata{i}, 1)/N;
    m(i,:) = mean(cdata{i});
  end

  v_w = 0; % variancia intraclasse (within class variance)
  for i = 1:k
    v_i = 0;
    for j = 1:size( cdata{i}, 1 )
      ssq = 0; % sum of square error
      for l = 1:d
        ssq = ssq + ( cdata{i}(j,l) - m(i,l) )^2; 
      end
      v_i = v_i + ssq;
    end
    v_w = v_w + w(1,i) * v_i;
  end

  m_t = 0; % media total
  for i = 1:N
    m_t = mean(data);
  end

  v_t = 0; % total variance
  
  for i = 1:N % total variance
    ssq = 0; % sum of square error
    for j = 1:d
      ssq = ssq + ( data(i,j) - m_t(j) )^2; 
    end
    v_t = v_t + ssq;
  end

  v_b = v_t - v_w;
  iccv = v_b/(v_b + v_w);

end