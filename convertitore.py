import Bio
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#Input e controlli sull'input
L1 = None
while L1 is None:       #L1 deve essere >= 1
    try:
        L1 = int(input("L1: "))
    except:
        print("Must be an integer")

L2 = None
while L2 is None or L2 <= L1:      #L2 deve essere maggiore di L1
    try:
        L2 = int(input("L2: "))
    except:
        print("Must be an integer")

Q1 = None
while Q1 is None:
    try:
        Q1 = float(input("Q1: "))
    except:
        print("Must be a number (float or int)")

Q2 = None
while Q2 is None or Q2 <= Q1:     #Q2 deve essere maggiore di Q1
    try:
        Q2 = float(input("Q2: "))
    except:
        print("Must be a number (float or int)")

P = -1
while P < 0 or P > 1:   #P deve essere compreso tra 0 e 1 (inclusi)
    try:
        P = float(input("P (in range[0, 1]): "))
    except:
        print("Must be a number in range [0, 1]")

file_input_name = input("File to convert: ")
fastq_records = SeqIO.parse(file_input_name, 'fastq')
fastq_record_list = list(fastq_records)
out = open('./output.fa', 'w')

#scansione dei read in input
for fastq_record in fastq_record_list:

    record_len = len(fastq_record.seq)
    phred_quality = fastq_record.letter_annotations['phred_quality']
    min_quality = min(phred_quality)
    
    #se il record fastq rispetta le richieste sulla lunghezza e sulla qualità minima, allora si può proseguire
    #altrimenti si passa al record successivo
    if record_len >= L1 and min_quality > Q1 and record_len <= L2:

        #calcolo delle sottoregioni regioni in cui ogni elemento della sottoregione ha qualità >= Q2
        bool_list = [quality >= Q2 for quality in phred_quality] #bool_list[i] = true se l'elemento i della sequenza ha qualità >= Q2
        #start_list[i] contiene l'indice di inizio dell'i-mo intervallo e end_list[i] contiene l'indice di inizio dell'i-mo intervallo
        start_list = [i for i in range(len(bool_list)) if bool_list[i] and (i == 0 or not bool_list[i - 1])] #contiene la posizione di inizio dei vari intervalli
        end_list = [i for i in range(len(bool_list)) if bool_list[i] and (i == len(bool_list) - 1 or not bool_list[i + 1])] #contiene la posizione di fine dei vari itervalli
        
        #si genera intervallo fittizio che inizia in posizione 2 e finisce in posizione 0
        #in questo modo se non esiste una sottoregione che ha basi con qualità maggiore o uguale a Q2,
        #il read viene scartato
        start_list[:0] = [2] 
        end_list[:0] = [0]

        #si calcola la lunghezza degli intervalli trovati e si sceglie l'intervallo con lunghezza + lunga
        intervals_length = [end_list[i] - start_list[i] + 1 for i in range(len(start_list))]
        max_breadth = max(intervals_length)
        index = intervals_length.index(max_breadth)
        start = start_list[index]
        end = end_list[index]

        #se la sottoregione è lunga almeno P% della lungheza del read allora si converte in fasta e si stampa in output
        #altrimenti si passa al read successivo
        #si passa al read successivo anche se non esiste una sottoregione, infatti in questo caso 
        #max_breadth = 0 - 2 + 1 = -1 e P * record_len è sempre maggiore o uguale a zero
        if max_breadth >= P * record_len: 
            #inserisce in description le informazioni necessarie
            description = "seq_len:" + str(record_len) 
            description += " min_quality:" + str(min_quality)
            description += " subregion_start_end:" + "[" + str(start) + ":" + str(end) + "]"

            mean = sum(phred_quality[start:end + 1])/max_breadth
            description += " subregion_medium_quality:" + str(mean)

            fastq_record.description = description

            fasta_record = fastq_record.format('fasta')

            #stampa su file di output il read convertito
            out.write(fasta_record)

out.close()