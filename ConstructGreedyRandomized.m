function x = ConstructGreedyRandomized(x, n, h, l, u, alpha, k, data, qtd, TimeToRun)

   S = 1:n;
   mem = ones(1,n);
   tam_S = n;
   epsilon = 0.0000001;
   
   while tam_S > 0 
      
      f_x = f(x, n, k, data, qtd);
      f_min = inf; f_max = -inf;
      z = zeros(1,n);
      
      for i = 1 : n
         if mem(i) == 1
            z(i,:) = LineSearch(x, h, i, n, l, u, k, data, qtd, f_x);
            
            g(i) = f(z(i,:), n, k, data, qtd);
            
            if f_min > g(i) 
               f_min = g(i); 
            end
            if f_max < g(i) 
               f_max = g(i);
            end
         end
         
      end
      
      RCL = [];
      for i = 1 : n
         if mem(i) == 1
            if(g(i) <= ((1 - alpha)*f_min + alpha*f_max) + epsilon)
               index = size(RCL);
               RCL(index(2)+1) = i;
            end
         end
      
      end
      
      tam = size(RCL);
      j = randi([1 tam(2)]);
      j = RCL(j);
      x(j) = z(j,j); 
      S(find(S == j)) = [];
      mem(j) = 0; 
      tam_S = tam_S - 1;
      
      if ttest(TimeToRun); return; end; % fim da execução?
      
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