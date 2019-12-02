# Kriging Ordinario

- Scrivere una funzione di octave che calcoli il il kriging ordinario in un insieme di punti P con coordinate note (Px, Py).

- Il prototipo della funzione sarà simile (ma non necessariamente uguale) al seguente:

`[pred, var] = ordkrig(D, P, semivarmod, par)`

dove:
* D = L’insieme dei valori noti in forma di matrice D [N x 3] dove N = numero dei dati. Per coerenza assumiamo che la 1a colonna di D rappresenti la coordinata x, la 2a colonna la coordinata y e la 3a colonna il valore z da interpolare.
* P = L’insieme degli M punti per i quali calcolare il kriging come matrice P [2 x M] con le coordinate x nella 1a colonna e y nella 2a.
* semivarmod =  modello di semivariogramma e i relativi parametri saranno passati con l’array “par”. Per generalizzare il più possibile usa un puntatore ad una funzione semivariogramma esterna. Non è necessario scrivere tutte i possibili modelli di semivariogramma, per provare la funzione sarà sufficiente usare uno dei modelli visti a lezione.

`par = parametri di semivarmod (nugget, range, sill, ecc)`

La funzione deve restituire un array  “pred” con il valore interpolato nei punti P e un altro array con la loro varianza “var”.

- Sperimentare la funzione ottenuta usando un set di dati a piacere fare un grafico dei dat interpolati e commentare i risultati ottenuti.

## Descrizione
La funzione ordkrig richiede in input i parametri nugget, sill e range. Quindi, per effettuare un test migliore, ho deciso di ricavarli dal semivariogramma e successivamente usarli nel test della funzione.

#### Operazioni preliminari
1. Utilizzo 100 campioni generati casualmente, contenuti in un file (campioni.dat) costituito da 3 colonne formattate come segue:

`Coordinata X – Coordinata Y – Misura`

2. Genero il semivariogramma per mezzo della funzione Octave “semivariogramma.m”.
3. Analizzo il semivariogramma regolarizzato con un modello parametrico sferico avente i seguenti valori:

- Nugget = 6000
- Sill = 2000
- Range = 150

![Grafico1](https://github.com/elpanas/Kriging-ordinary/blob/master/Grafici/2.png)

Effettuo il fitting pesato e verifico che non è necessario effettuare un detrend (μ(x) costante)
Proseguo quindi con il kriging ordinario.

### Kriging ordinario
Mettiamo i 100 campioni casuali D su un grafico per capirne meglio la posizione.
Scegliamo dei punti P di cui non conosciamo la misura associata. Da sinistra a destra abbiamo:

 P1 | P2 | P3 
| -------- | -------- | --------- |
| (170;40) | (230;50) | (390;100) |

Utilizziamo la funzione ordkrig() per calcolare i previsori e la loro rispettiva varianza.

 P(x;y)  | Previsore P* | Varianza sigma2 
|---------|:---------------:|:-----:|
| P1(170;40) | 902.26 | 6000.7 |
| P2(230;50) | 935.01 | 6002.1 |
| P3(390;100) | 981.58 | 6005.4 |

Metto i risultati sul grafico:

![Grafico2](https://github.com/elpanas/Kriging-ordinary/blob/master/Grafici/3.png)

### Conclusioni
Possiamo notare come la varianza del valore al punto P3 sia superiore a quella degli altri punti, a causa della sua posizione maggiormente isolata; lo stesso comportamento del valore in posizione P2. Al contrario il valore nel punto P1 ha la varianza più bassa dei 3.
Abbiamo quindi dimostrato come il previsore in punti vicini o in gruppo sia più preciso (con una varianza più bassa) rispetto alla previsione in punti più isolati.
