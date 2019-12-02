function [pred,var] = ordkrig(D, P, model, par)
  %{
	Input:      D - coordinate e valori dei punti noti
              P - coordinate dei punti con valori ignoti
          model - modello di semivariogramma
            par - array contenete valori di semivarianza
                   sill - soglia oltre la quale la semivarianza non aumenta
                  range - valore oltre il quale la varianza non aumenta
                          o distanza max di P. Cioè oltre questo limite non
                          ci sono valori di P.
                 nugget - è l'errore di osservazione in misure ripetute
                          o un margine nel semivariogramma (effetto "pepita")

    Output: pred - vettore contenente i valori dei punti ignoti (previsori)
             var - vettore contenente le varianze dei previsori

	Obiettivo: la funzione ordkrig esegue il kriging ordinario, prendendo
             in input i punti noti (coordinate e valore) e restituendo
             in output il valore stimato o previsore (P_star) nei punti P
             e la rispettiva varianza nei punti.

  --- TEORIA ALLA BASE DEL KRIGING ---
  Modello geostatistico: G(x) = u(x) + e(x)
	definito da:
		u(x) - componente su larga scala (trend)
		e(x) - componente su piccola scala (stocastica stazionaria)
	Di questo modello fanno parte sia i punti noti D, che quelli ignoti P.
	Studiando come variano i valori tra i campioni (D) posti ad una data
	distanza da P, è possibile prevedere o stimare il valore in un punto vicino (P).
	La variazione suddetta è detta semivarianza, espressa in funzione della
	distanza dei punti per i quali è calcolata.
	Il grafico della funzione semivarianza è il semivariogramma.

	Modello parametrico: è una funzione che serve a regolarizzare il semivariogramma
						 empirico per avere
						 - risultati più affidabili
						 - varianza non negativa
	In questo caso è stato scelto il modello sferico.
  ------------------------

  --- KRIGING ---  
	Per stimare il valore nel punto/i P esistono diverse tecniche, tra cui il Kriging.

	Kriging: il kriging è una tecnica geostatistica con cui si stima
			     un valore/misura (Previsore) in un punto P di coordinate (Px,Py),
           sulla base dei valori noti di altri punti D, posti ad una data
           distanza da quello cercato, che viene interpolato a partire
           dai valori noti.
           Per stimare il valore nel punto/i di interesse P usiamo i punti D adiacenti.
           a cui assegnamo dei pesi, in modo da dare a questi ultimi più o meno importanza.
           Usiamo però un approccio geostatistico. Cioè i pesi vengono ricavati dalle
           proprietà statistiche dei punti noti che si assumono essere stazionari.
           Usa il metodo dei pesi per dare maggiore importanza (peso) ai punti
           vicini a quello cercato, rispetto a quelli lontani o più isolati.
           Inoltre fornisce la varianza del valore stimato, cioè l'errore o per
           meglio dire: quanto si discosta dal valor medio u(x).
           Anche per questo motivo è importante che u(x) sia costante.

	K. ordinario: è la versione del kriging nel caso di u(x) costante.
                Il previsore è dato dalla moltiplicazione dei pesi per i
                valori dei rispettivi punti. Posto che la somma dei pesi
                deve essere uguale a 1, i pesi migliori (cioè quelli che
                minimizzano gli errori) sono calcolati con il metodo del
                moltiplicatore di Lagrange.
                Ricaviamo i valori che ci servono dalla seguente formula.
                La formula è: 		  gammaDD * w + lambda = gammaPD
                da cui si ricava: 	w + lambda = gammaDD / gammaPD

			dove:	      w - vettore dei pesi (con lambda in ultima posizione)
             lambda - moltiplicatore di Lagrange
            gammaDD - matrice delle semivarianze tra i punti noti D
            gammaPD - matrice delle semivarianze tra i punti P e quelli
                      noti D (con il coeff. "1")
                      
            Semplificando: se esistono proprietà statistiche simili tra i punti
            noti D, essendo i punti P vicini, c'è un'alta probabilità che le
            stesse proprietà sussistano tra i punti P e quelli noti D. Sulla base 
            di questa somiglianza, possiamo stimare D.

            Calcola il previsore (P_star) cioè la stima del valore nei punti P
            La formula è: P_star = w' * V
            dove: w'- vettore trasposto dei pesi
                  V - vettore con i valori dei punti noti a cui associare i pesi

            Infine calcola la varianza (sigma2) dei punti P, definita come
            la varianza della differenza del previsore pesato G_star meno
            il valore in P.
            La formula è: sigma2 = lambda + w' * gammaPD_0
            dove: 	sigma2 - varianza nel punto P (o errore previsto)
                    lambda - moltiplicatore di Lagrange
                        w' - vettore trasposto dei pesi
                 gammaPD_0 - matrice delle semivarianze tra i punti ignoti P
                             e i punti noti D (senza il coeff. "1")
	-------------------------------------------------------------------------
  %}

  % --- INIZIALIZZA LE VARIABILI ---
  nugget = par(1);    % errore di misura (effetto "pepita")
  sill = par(2);      % soglia massima della semivarianza
  range = par(3);     % valore oltre il quale la varianza non aumenta o distanza max di P
  cp = 0;             % contatore dei valori ignoti dei punti P
  i = 1;              % indice dei vettori di output
  [rp,Mp] = size(P);  % numero di punti ignoti P (righe di P)
  % --------------------------------

  % --- SEMIVARIANZA TRA I PUNTI NOTI D ---
  %{
	ogni riga contiene una tripla di valori (le coordinate Dx,Dy e la misura)
	che corrispondono ad un punto. Quindi il numero di righe è anche il numero
	dei punti campionati.
  %}
  [rd,Nd] = size(D); % numero di punti con coordinate note (righe "rd")

  % scorro le righe e quindi i punti noti
  for n = 1:rd % n indice al punto noto attuale
    delDD(:,n,:) = abs(D(n,:)-D); % matrice con le differenze
  endfor

  % calcola le distanze euclidee tra i punti noti (lag tra D1 e D2)
  lagDD = sqrt(delDD(:,:,1).^2 + delDD(:,:,2).^2);

  %{
	calcola la semivarianza delle distanze tra i punti noti
	sulla base dei valori di sill e range forniti in input
	con modello parametrico di tipo Sferico.
  %}
  switch (model)
    case 'gaussian' % puntatore a funzione
      gammaDD = semivarfunc(lagDD, par, @gaussian_model); 
    case 'spherical' % puntatore a funzione
      gammaDD = semivarfunc(lagDD, par, @spherical_model); 
  endswitch

  V = D(:,3); % vettore dei valori nei punti noti

  % --- COEFFICIENTI PER IL MOLTIPLICATORE DI LAGRANGE... ---

  % aggiungo i coeff. alla semivarianza tra i punti noti D
  %{
	E' necessario aggiungere una riga e una colonna di "1" per il
	moltiplicatore di Lagrange, che sarà nel vettore con i pesi.
	Ad esempio: se abbiamo 5 campioni D, il vettore dei pesi sarà 
	lungo 5+1 che è il posto del moltiplicatore di Lagrange; quindi
	per eseguire l'operazione successiva, le matrici devono avere
	uguale numero di righe e colonne.
  %}
  gammaDD(:,end+1) = 1; % colonna di "1"
  gammaDD(end+1,:) = 1; % riga di "1"
  gammaDD(end,end) = 0; % valore "0"
  % -----------------------

  % --- CICLA I PUNTI 'P' E NE CALCOLA IL VALORE E LA VARIANZA ---

  for cp = 1:rp % ciclo i punti con valore sconosciuto
      % --- SEMIVARIANZA TRA IL PUNTO P IGNOTO E I PUNTI D NOTI ---
      %{
      la distanza euclidea di punti con 2 coordinate si calcola con la radice di:
      (la differenza delle coordinate x) + (quella delle coordinate y), al quadrato
      %}
      % uso la funzione abs nel caso in cui la coordinata x1 < x2 (opp. y1 < y2)
      delPD = abs(P(cp,:)-D(:,[1,2]));
      lagPD = sqrt(delPD(:,1).^2 + delPD(:,2).^2); % distanza tra P e D

      %{
      calcola la semivarianza delle distanze PD sulla base
      dei valori di sill e range forniti in input, con modello
      parametrico di tipo Sferico.
      %}
      switch (model)
        case 'gaussian'
          gammaPD_0 = semivarfunc(lagPD, par, @gaussian_model);
        case 'spherical'
          gammaPD_0 = semivarfunc(lagPD, par, @spherical_model);
      endswitch

      % --- COEFFICIENTI PER IL MOLTIPLICATORE DI LAGRANGE... ---
      %{
      aggiungo il coeff. alla matrice delle semivarianze tra il punto P
      e i punti D noti.
      Salvo il vettore originale per evitare di modificarlo, perchè
      servirà  nel calcolo dell'errore sigma2.
      %}
      gammaPD = gammaPD_0;
      gammaPD(end+1) = 1; % matrice con l'aggiunta del valore "1"

      % --- PESI ---
      w = gammaDD\gammaPD; % calcola la matrice con i pesi

      % --- PREVISORE P* ---
      P_star = w(1:rd)' * V; % stima del punto

      % --- VARIANZA NEL PUNTO ---
      %{
       w - vettore dei pesi
      rd - numero di righe, cioè di punti D noti sulla base dei quali P è stimato
      %}
      lambda = w(rd+1); % estrae il moltiplicatore di Lagrange dall'ultima posizione del vettore w
      sigma2 = lambda + w(1:rd)' * gammaPD_0; % vettore delle varianze dei rispettivi punti P

      % --- MEMORIZZA I RISULTATI NEI VETTORI ---
      pred(i,:) = P_star; % vettore dei previsori o valori interpolati
      var(i,:) = sigma2; % vettore delle varianze dei previsori
      % -----------------------

      i = i + 1; % incrementa l'indice dei vettori di output
  endfor

endfunction
