# Kriging Ordinario
Progetto di esame Luca Panariello 27/05/2019


- Scrivere una funzione di octave che calcoli il il kriging ordinario in un insieme di punti P con coordinate note (Px, Py).

- Il prototipo della funzione sarà simile (ma non necessariamente uguale) al seguente:

[pred, var] = ordkrig(D, P, semivarmod, par)

dove:
D = L’insieme dei valori noti in forma di matrice D [N x 3] dove N = numero dei dati. Per coerenza assumiamo che la 1a colonna di D rappresenti la coordinata x, la 2a colonna la coordinata y e la 3a colonna il valore z da interpolare.
P = L’insieme degli M punti per i quali calcolare il kriging come matrice P [2 x M] con le coordinate x nella 1a colonna e y nella 2a.
semivarmod =  modello di semivariogramma e i relativi parametri saranno passati con l’array “par”. Per generalizzare il più possibile usa un puntatore ad una funzione semivariogramma esterna. Non è necessario scrivere tutte i possibili modelli di semivariogramma, per provare la funzione sarà sufficiente usare uno dei modelli visti a lezione.
par = parametri di semivarmod (nugget, range, sill, ecc).

La funzione deve restituire un array  “pred” con il valore interpolato nei punti P e un altro array con la loro varianza “var”.

- Sperimentare la funzione ottenuta usando un set di dati a piacere fare un grafico dei dat interpolati e commentare i risultati ottenuti.
