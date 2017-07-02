function [x_best fconv nimprov] = LocalSearch(x, n, h, l, u, MaxDirToTry, k, data, qtd, fconv, nimprov, TimeToRun)

   Improved = true; D = [];
   x_best = x; f_best = f(x, n, k, data, qtd);
   NumDirToTry = min(power(3, n) - 1, MaxDirToTry);
   
   ctime = toc;
   firstImp = true;

   while Improved
   
      Improved = false;
      
      while (length(D) <= NumDirToTry) && (~Improved) 
         
         if length(D) == NumDirToTry
            break;
         end
         
         r = 0;
         while ismember(r, D) || r == 0
            if n > 32 % limite para a funcao randi
               r = randi([1 (power(3, 32) - 1)]);
            else
               r = randi([1 (power(3, n) - 1)]);
            end;
         end
         
         index = size(D);
         D(index(2)+1) = r;
         
         d = Ternary(r, n);
         
         x = x_best + h*d; % calcula novo ponto dado o vetor direção e a densidade da grade de busca h

         if (l <= x) 
            if (x <= u)
               fo = f(x, n, k, data, qtd);
               if fo < f_best
                  x_best = x; f_best = fo;
                  if toc - ctime > 5 || firstImp
                    fconv(nimprov, 1) = toc;
                    fconv(nimprov, 2) = f_best;
                    nimprov = nimprov + 1;
                    ctime = toc;
                    firstImp = false;
                  end
                  D = [];
                  Improved = true;
               end
            end
         end
      
         if ttest(TimeToRun); return; end; % fim da execução?
      
      end

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