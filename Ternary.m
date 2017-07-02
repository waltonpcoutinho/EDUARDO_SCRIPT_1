function d = Ternary(r, n)

   str = dec2base(r, 3, n);
   str = str2double(regexp(num2str(str), '\d','match'));
   d = zeros(1, n);
   
   for i=1:n
      if str(i) == 2
         d(i) = -1;
      else
         d(i) = str(i);
      end
   end

end
