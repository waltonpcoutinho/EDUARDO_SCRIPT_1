function [x ImprC] = ConstructGreedyRandomizedImproved(x, n, h, l, u, k, data, qtd, TimeToRun)

   ImprC = false;

   S = 1:n;
   alpha = rand();
   ReUse = false;   
   mem = ones(1,n);
   tam_S = n;
   epsilon = 0.0000001;
   
   z = zeros(1,n); 
   f_min = inf; f_max = -inf;
   f_x = f(x, n, k, data, qtd);
   
   while tam_S > 0 
      
      for i = 1 : n
         if mem(i) == 1

            if ReUse == false
                z(i,:) = LineSearch(x, h, i, n, l, u, k, data, qtd, f_x);
                g(i) = f(z(i,:), n, k, data, qtd);
            end
            
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
      
      if x(j) == z(j,j) 
        ReUse = true;
      else
        x(j) = z(j,j); 
        ReUse = false;
        if g(j) < f_x
          f_x = g(j);
          ImprC = true; 
        end
        f_min = inf; f_max = -inf;
      end
      
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