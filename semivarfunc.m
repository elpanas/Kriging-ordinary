function gamma = semivarfunc(lag, par, modelf)
  
  %{
    Input:     lag - distanza euclidea tra i punti
               par - parametri per il modello di semivariogramma
            modelf - funzione modello di semivariogramma
    
    Output:  gamma - semivarianza;
  %}
  
  nugget = par(1);  % errore di misura ("pepita")
  sill = par(2);    % soglia massima di semivarianza
  range = par(3);   % valore oltre il quale la varianza non aumenta
  
  % calcolo della semivarianza
  gamma = nugget + modelf(lag, [sill,range]);
  
endfunction
