function Bh_x = RandomlySelectElement(x,n,h,l,u)

    %{ 
    if n > 32 % limite para a funcao randi (DIREÇÕES POSSÍVEIS INCOMPLETAS)
      r = randi([1 (power(3, 32) - 1)]);
    else
      r = randi([1 (power(3, n) - 1)]);
    end;
    %}
    
    % Escolhe aleatoriamente um ponto do Grid de tamanho de passo h.
    for i=1:n
	    % Calcula o numero de pontos naquela direcao com tamanho de passo h.
	    numPointsPos = floor((u(i) - x(i))/h);
	    numPointsNeg = floor((x(i) - l(i))/h);
	    numPoints = numPointsPos + numPointsNeg;
	    %(u(i) - x(i))/h
        if numPoints == 0
            xGrid(i) = x(i); 
        else
            % Escolhe aleatoriamente um dos indices do Grid na direcao i.
            gridI = (randi([0 numPoints]) - numPointsNeg); 
            % Calcula a posicao da direcao no eixo i.
            xGrid(i) = x(i) + gridI*h;
        end
    end

    % Calcula a projecao do vetor de x a xGrid.
    normaVector = norm(x - xGrid); 
    if normaVector ~= 0
       for i=1:n
           Bh_x(i) = x(i) + ((h* (xGrid(i)- x(i)))/normaVector);
       end
    else
        Bh_x = x;
    end
   
    %{
    d = Ternary(r, n);
    d = d*h;
    Bh_x = x + d;
    Bh_x = x + h*(Bh_x - x)/norm(Bh_x - x);
    %}
   
end
