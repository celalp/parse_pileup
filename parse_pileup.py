import pandas as pd
import argparse
import re
from collections import Counter
import logging as log


#############################################
# this is script processes a pileup file and
# returns the best call for each nucleotide
# this is designed for tiny "genomes" such 
# as the mitochondiral genome not intended
# for larger variant calls 
############################################

def variant_counter(pileup, reference):
    """
    This function is responsible for getting the variant calls for a given nucleotide
    :param pileup: The pileup column of the file this is passed within the read_pileup function
    :param reference: iterable reference sequence
    :return: a dataframe of counts, this is merged inside the read_pileup function
    """
    counts = {"A": [], "T": [], "G": [], "C": [], "*":[], "ins": [],
              "del": [], "ins_type": [], "del_type": []}
    for k in range(len(pileup)):
        ref=reference[k]
        variants=pileup[k]
        if len(variants)==0:
            for key in counts.keys():
                counts[key].append(0)
        else:
            search = {"A": 0, "T": 0, "G": 0, "C": 0, "*":0, "ins": 0,
                      "del": 0, "ins_type": [], "del_type": []}
            i=0
            while i < len(variants):
                j = 1
                if variants[i] in ".,":
                    search[ref]+=1
                    i += 1
                elif variants[i] in "Aa":
                    search["A"] +=1
                    i += 1
                elif variants[i] in "Gg":
                    search["G"] +=1
                    i += 1
                elif variants[i] in "Cc":
                    search["C"] +=1
                    i += 1
                elif variants[i] in "Tt":
                    search["T"] +=1
                    i += 1
                elif variants[i] is "*":
                    search["*"]+=1
                    i+=1
                elif variants[i] is "-" and variants[i+j] in "0123456789" and i+j+1<len(variants): # deletions
                    # get the number of digits
                    while variants[i+j] in "0123456789":
                        j+=1
                    deletion_length=int(variants[i+1:i+j]) #convert to int to get the actual number
                    search["del"]+=1
                    search["del_type"].append(variants[i+j:i+j+deletion_length])
                    i+=1+j+deletion_length
                elif variants[i] is "+" and variants[i+j] in "0123456789" and i+j+1<=len(variants): # insertions same as deletions could be dryer
                    while variants[i+j] in "0123456789":
                        j+=1
                    insertion_length=int(variants[i+1:i+j])
                    search["ins"]+=1
                    search["ins_type"].append(variants[i+j : i+j+insertion_length])
                    i+=1+j+insertion_length
                else:
                    i+=1
            for key in search.keys():
                if key is not "ins_type" and key is not "del_type":
                    counts[key].append(search[key])
                elif key is "ins_type" or key is "del_type":
                    if len(search[key])>0:
                        upper=[x.upper() for x in search[key]]
                        to_append=Counter(upper).most_common(1)[0][0]
                        counts[key].append(to_append)
                    else:
                        counts[key].append(None)
    counts=pd.DataFrame(counts)
    return counts

def generate_bestcall(df):
    """
    This function takes the results of the parse_pileup function and figures out what the best call is
    if the best call is an insertion then the next best call is used. What this function DOES NOT do
    is to check if the next best call has another call with the same number of reads. If there are no
    reads mapping to the function then _ is used to indicate there is nothing there
    :param df: results of parse_pileup
    :return: a list of the best calls
    """
    bestcall=df.iloc[:,:6].idxmax(axis=1)
    not_ins=df.iloc[:,:5].idxmax(axis=1)
    actual_call=[]
    for i in range(len(bestcall)):
        if df.iloc[i,:][bestcall[i]]>0 :
            if bestcall[i] is not "ins" and bestcall[i] is not "*":
                actual_call.append(bestcall[i])
            elif bestcall[i] is "*":
                actual_call.append("")
            elif bestcall[i] is "ins": #now i need to append this to the next abundant call
                next_best=not_ins[i]
                if next_best is "*":
                    next_best=""
                if df.iloc[i,:][next_best]>0:
                    actual_call.append(next_best+bestcall[i])
                else:
                    actual_call.append("_")
                    msg="The best match was an insertion but there are 0 reads as second best match at position "+str(i)
                    log.warning(msg)
        else:
            actual_call.append("_")
            msg="There are 0 reads at position "+str(i)
            log.warning(msg)
    return actual_call

def heteroplasmy(A, T, G, C, deletion):
    """
    heteroplasmy calculation, only uses A,T,G,C and deletion counts
    :param A: A count
    :param T: T count
    :param G: G count
    :param C: C coung
    :param deletion: * count
    :return: a list of float
    """
    if sum([A, T, G, C, deletion])==0:
        het=None
    else:
        het=1-(max([A, T, G, C, deletion])/sum([A, T, G, C, deletion]))
    return het

def read_pileup(filename, ref_fasta):
    """
    This is the workhorse function that calls all the other functions above
    :param filename: pileup file
    :param ref_fasta: reference fasta
    :return: a dataframe that is to be written to file
    """
    pileup_info=pd.read_csv(filename, sep="\t", header=None, usecols=[0,1,2,3,4],
                            names=["Ref", "Pos", "RefCall", "Depth", "Pileup"])
    pileup_info["Ref"] = "rCRS"
    for i in range(pileup_info.shape[0]):
        if pileup_info.iloc[i,1]>144:
            pileup_info.iloc[i, 1]=pileup_info.iloc[i,1]-144
        else:
            pileup_info.iloc[i, 1] = pileup_info.iloc[i, 1] + 16424
    ref=[]
    for line in open(ref_fasta):
        if ">" in line:
            continue
        else:
            line = line.rstrip("\n")
            ref.append(line)

    reference="".join(ref)
    reference=reference[144:]+reference[:144]
    merger=pd.DataFrame({"Ref":["rCRS"]*16568, "Pos":list(range(1,16569))})
    pileup_info=merger.merge(pileup_info, how="outer", on=["Ref", "Pos"])
    pileup_info["RefCall"]=list(reference)
    pileup_info["Depth"].fillna(0, inplace=True)
    pileup_info["Pileup"].fillna(value="", inplace=True)
    variant_counts=variant_counter(pileup_info["Pileup"], reference)
    pileup_info=pileup_info.drop(["Pileup"], axis=1).join(pd.DataFrame(variant_counts))
    pileup_info["BestCall"]=generate_bestcall(variant_counts)
    pileup_info["Heteroplasmy"]=pileup_info.apply(lambda x: heteroplasmy(x['A'], x['T'], x['G'], x['C'], x['*']), axis=1)
    return pileup_info

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='parse an mpileup file')
    parser.add_argument('-m', '--mpileup', type=str, help='mpileup file', action="store")
    parser.add_argument('-r', '--reference', type=str, help="reference fasta")
    parser.add_argument('-f', '--fasta', help="genererate separate fasta file instead of appending to the end of the file",
                        action='store_true')
    parser.add_argument('-o', '--output', type=str, help='name of the output file', action="store")
    args = parser.parse_args()

    pileup=read_pileup(args.mpileup, args.reference)
    pileup.to_csv(args.output, index=False)
    fasta = "".join(pileup["BestCall"])
    if args.fasta:
        fastaname=args.output+".fa"
        with open(fastaname, "w") as f:
            f.write((">"+pileup["Ref"][0]+"\n"))
            f.write(fasta)
            f.close()
    else:
        with open(args.output, "a") as f:
            f.write(fasta)
            f.close()





