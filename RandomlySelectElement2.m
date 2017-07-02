function Bh_x = RandomlySelectElement2(x,n,h,l,u)

    if n > 32 % limite para a funcao randi (DIREÇÕES POSSÍVEIS INCOMPLETAS)
      r = randi([1 (power(3, 32) - 1)]);
    else
      r = randi([1 (power(3, n) - 1)]);
    end;
    
    d = Ternary(r, n);
    d = d*h;
    Bh_x = x + d;
    Bh_x = x + h*(Bh_x - x)/norm(Bh_x - x);
   
end
