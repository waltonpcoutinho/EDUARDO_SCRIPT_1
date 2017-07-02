%{
function value = f(x, n, k, data, qtd) 

   soma_total = 0; 
   D = n/k;
   
   for o=1:qtd;
      dim = 0;
      mi = inf;
      for c = 1:k;
          soma_c = 0;
          for j = 1:D;
              dim = dim + 1;
              soma_c = soma_c + (data(o,j) - x(dim))^2;
          end;    
          dist_c(c) = sqrt(soma_c);
          if dist_c(c) < mi
            mi = dist_c(c);
          end;
      end; 
      soma_total = soma_total + mi;
   end;
   
   value = soma_total;
   
end
%}
function value = f(x, n, k, data, qtd) 

   soma_total = 0; 
   D = n/k;
   
   for o=1:qtd;
      dim = 0;
      mi = inf;
      for c = 1:k;
          soma_c = 0;
          for j = 1:D;
              dim = dim + 1;
              soma_c = soma_c + (data(o,j) - x(dim))^2;
          end;    
          dist_c(c) = sqrt(soma_c);
          if dist_c(c) < mi
            mi = dist_c(c);
          end;
      end; 
      soma_total = soma_total + mi;
   end;
   
   value = soma_total;
   
end
