# Convertitore da Fastq a Fasta
## Progetto Elementi di Bioinformatica - UniMiB

### Input
Il convertitore richiede in input da console i seguenti parametri:
  * L1 -> Intero che rappresenta la lunghezza minima del read
  * L2 -> Intero (maggiore di L1) che rappresenta la lunghezza massima del read
  * Q1 -> Numero che rappresenta la qualità minima delle basi che compongono il read
  * Q2 -> Numero (maggiore di Q2) che rappresenta la qualità minima richiesta di una sottoregione della sequenza lunga almeno P% della lunghezza totale del read
  * P -> Numero compreso tra 0 e 1 (inclusi) che rappresenta la percentuale della lunghezza minima consentita per la sottoregione le cui basi hanno qualità maggiore o uguale a Q2
  * input_file -> Stringa che contiene il nome del file da cui il convertitore leggerà i read in formato fastq
  
### Output
Il convertitore stamperà in output sul file "output.fa" i read tradotti in formato fasta che rispettano i vincoli, nella descrizione del read vengono inserite le seguenti informazioni:
  * Lunghezza del read tradotto
  * Qualità minima del read
  * L'indice di inizio e l'indice di fine della sottoregione le cui basi hano qualità maggiore o uguale a Q2
  * La qualità media della sottoregione
  
### Note
  * Tra tutte le sottoregioni valide si sceglie quella con ampiezza massima
  * Se non esiste una sottoregione del read che ha basi con qualità maggiore o uguale a Q2, allora il read non viene stampato in output
