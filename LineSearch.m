function z = LineSearch(x, h, i, n, l, u, k, data, qtd, f_x)

   z = x;
   
   z1 = x; z2 = x;
   v = x(i);
   z1(i) = v + h;
   z2(i) = v - h;
   
   best_f = f_x;
   
   if (z1 >= l) & (z1 <= u)
      f_z1 = f(z1, n, k, data, qtd);
      if f_z1 < best_f
         z = z1;
         best_f = f_z1;
      end
   end

   if (z2 >= l) & (z2 <= u)
      f_z2 = f(z2, n, k, data, qtd);
      if f_z2 < best_f
         z = z2;
      end
   end
   
end
