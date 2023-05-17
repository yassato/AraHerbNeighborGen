########################################################################
# Python script to compile A. lyrata and A. thaliana alleles from .maf #
########################################################################

# load library and data
import pandas as pd
import re

def spl(word):
    return [char for char in word]

file1 = open('arabidopsis_thaliana_TAIR10_vs_arabidopsis_lyrata_v_1_0_lastz_net.chr5.maf', 'r') # change the input file when compiling another chr.
Lines = file1.readlines()

# search position and sequence
pos_seq = []
Atha_seq = []
Alyr_seq = []
score_seq = []
for i in range(0,len(Lines)):
    MatchLine = re.search("arabidopsis_thaliana",Lines[i])
    if MatchLine != None:
        print(i)
        print(len(Lines[i].split()[6]) == len(Lines[i+1].split()[6]))

        #Atha
        tha = Lines[i].split()
        pos = int(tha[2])
        seq = tha[6]
        splt = spl(tha[6])

        seqLen = len(splt)
        start = pos
        end = pos +seqLen

        for k in splt:
            Atha_seq.append(k)


        for j in range(start,end):
            pos_seq.append(j)

        #Alyr
        lyr = Lines[i+1].split()
        seq2 = lyr[6]
        splt2 = spl(lyr[6])

        for l in splt2:
            Alyr_seq.append(l)

        # score
        srch = re.search("[0-9]*\n",Lines[i-1])
        score = int(str(Lines[i-1])[srch.start():(srch.end()-1)])
        for m in range(seqLen):
            score_seq.append(score)


# write results
out = pd.DataFrame(index=[],columns=[])
out["Pos"] = pos_seq
out["Alyr"] = Alyr_seq
out["Atha"] = Atha_seq
out["Score"] = score_seq
out = out[out.Atha != out.Alyr]
out.to_csv("chr5.csv") # change the name of output when targeting another chr.

        
