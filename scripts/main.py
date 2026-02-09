import random
#simu53_BRCA#
human_p53_coding = open('data/human_p53_coding.txt','r')
codon_AA = open('data/codon_AA.txt','r')
#functions
def file_to_str(sequence):#מהקובץ לסטרינג
    DNA_strand=""
    for line in sequence:
        if line.startswith(">"):#מדלג על השורה הראשונה
            continue
        DNA_strand += line.strip()
    return DNA_strand

def DNA_RNA_Cod(DNA_strand):#פונקציה שמדמה את תהליך התעתוק ל- RNA
  RNA_strand=[]
  for line in DNA_strand:
    if line.startswith(">"):#מדלג על השורה הראשונה
      continue
    else:
      for i in line:#מעביר ל RNA
        i = i.upper()
        if i=="T":
            RNA_strand.append("U")
        else:
            RNA_strand.append(i) 
  return "".join(RNA_strand)


def Read_dict(codon_AA):#פונקציה שקוראת את המיפוי של הקודונים וחומצות האמינו מהקובץ codon_AA.txt ומציבה ב- dictionary.
  RNA_codon_table = {}#משתנה גלובלי
  for line in codon_AA:
    line = line.strip()
    if line == "":
        continue
    codon, amino_acid = line.split()
    RNA_codon_table[codon]=amino_acid
  return RNA_codon_table

def RNA_prot(RNA_strand, RNA_codon_table): #מתרגמת את רצף ה- RNA לרצף חלבון
    protein = []
    i = 0
    stop_codons = {"UGA", "UAG", "UAA"}

    while i + 3 <= len(RNA_strand):
        codon = RNA_strand[i:i+3]
        if codon in stop_codons:
            break
        if codon in RNA_codon_table:
            protein.append(RNA_codon_table[codon])
            i += 3
        else:
            i += 1

    return "".join(protein)

def Mutate_DNA(DNA_strand):#מוטצית החלפה
    DNA_strand = DNA_strand.upper()
    mutated_DNA = list(DNA_strand)#נעביר לליסט
    length= len(DNA_strand)
    placement = random.randrange(0, length)#מיקום המוטציה
    base = mutated_DNA[placement]
#------------------------------------------------       
    if base == "A":#כדי לוודא שזה משתנה והמוטציה לא משתנה לאותו האות
        random_TGC = random.randrange(1, 4)
        if random_TGC == 1:
            mutated_DNA[placement] = "T"
        elif random_TGC == 2:
            mutated_DNA[placement] = "C"
        elif random_TGC == 3:
            mutated_DNA[placement] = "G"
#------------------------------------------------         
    elif base == "T":
        random_AGC = random.randrange(1, 4)
        if random_AGC == 1:
            mutated_DNA[placement] = "A"
        elif random_AGC == 2:
            mutated_DNA[placement] = "G"
        elif random_AGC == 3:
            mutated_DNA[placement] = "C"
    #------------------------------------------------         
    elif base == "C":
        random_ATG = random.randrange(1, 4)
        if random_ATG == 1:
            mutated_DNA[placement] = "A"
        elif random_ATG == 2:
            mutated_DNA[placement] = "T"
        elif random_ATG == 3:
            mutated_DNA[placement] = "G"
    #------------------------------------------------         
    elif base == "G":
        random_ATC = random.randrange(1, 4)
        if random_ATC == 1:
            mutated_DNA[placement] = "A"
        elif random_ATC == 2:
            mutated_DNA[placement] = "T"
        elif random_ATC == 3:
            mutated_DNA[placement] = "C"
    mutated_DNA = "".join(mutated_DNA)
    return mutated_DNA

def insertion(DNA_strand):#מוטצית הוספה
    length= len(DNA_strand)
    placement = random.randrange(0, length+1)#מיקום המוטציה
    ACTG_base=""
    amount_to_add=random.randrange(1,4)
    counter=0
    while counter < amount_to_add:         
        ACTG_base_random = random.choice(["A", "T", "C", "G"])
        ACTG_base += ACTG_base_random
        counter += 1
    mutated_DNA = DNA_strand[:placement] + ACTG_base + DNA_strand[placement:]
    return mutated_DNA 

def deletion(DNA_strand):  # מוטצית החסרה
    length = len(DNA_strand)
    amount_to_delete = random.randrange(1, 4)  
    if length <= amount_to_delete:
        return DNA_strand
    placement = random.randrange(0, length - amount_to_delete + 1)  # מיקום המחיקה
    mutated_DNA = DNA_strand[:placement] + DNA_strand[placement + amount_to_delete:]
    return mutated_DNA

def Comp_seq(old,new):#השוואה בין החלבון הקודם לחדש אחרי מוטציה
  old_new_the_same=0
  old_new_the_diffrent=0
  i=0
  while i < len(old) and i < len(new):
    if old[i] != new[i]:
        old_new_the_diffrent += 1
    else:
        old_new_the_same += 1
    i += 1
  old_new_the_diffrent += abs(len(old) - len(new))
  return old_new_the_diffrent





#main program
Y_OR_N=input("Does the woman have a family genetic background and does she have one mutation in BRCA1,2 or not? (Y=Yes, N=No)")

if Y_OR_N == "Y":
    mutations_needed = 1
else:
    mutations_needed = 2

DNA_original = file_to_str(human_p53_coding)
#dictinary
RNA_codon_table = Read_dict(codon_AA)
#original protein
RNA_original = DNA_RNA_Cod(DNA_original)
protein_original = RNA_prot(RNA_original, RNA_codon_table)
#1000 times simulations
runs = 1000
generations_list = []
mutation_rate= 1e-4
for run in range(runs):
    gen = 0
    hits = 0

    DNA_current = DNA_original
    protein_current = protein_original

    while hits < mutations_needed:
        gen += 1
        if random.random() > mutation_rate:#אם אין מוטציות בדור הזה נמשיך להבא
            continue
        random_percent = random.randrange(1,101) 
        if random_percent <= 98:  # החלפת בסיס
            DNA_mut = Mutate_DNA(DNA_current)
        elif random_percent == 99:  # הוספה של בסיס
            DNA_mut = insertion(DNA_current)
        else:  # מחיקה של בסיס
            DNA_mut = deletion(DNA_current)
        #תרגום לחלבון
        RNA_mut = DNA_RNA_Cod(DNA_mut)
        protein_mut = RNA_prot(RNA_mut, RNA_codon_table)

        if Comp_seq(protein_current, protein_mut) > 0:#בודק אם החלבון השתנה ביחס לחלבון הנוכחי
            hits+=1
            protein_current = protein_mut
        
        DNA_current = DNA_mut#עובר לדור הבא עם הדנא החדש
    generations_list.append(gen)


#ממוצע והדפסה 
if Y_OR_N == "Y":
    print("For a female that does have BRCA1,2 Mutation:")
elif Y_OR_N == "N":
    print("For a female that does/ doesn't not have BRCA1,2 Mutation:")

avg_generations = sum(generations_list) / len(generations_list)
avg_years = avg_generations / 365
print("The mutation that will change the P53 protein will take in average %.2f years." %avg_years)

#סגירת קצבים (:
human_p53_coding.close()
codon_AA.close()