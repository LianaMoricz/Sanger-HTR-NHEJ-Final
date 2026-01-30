from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import sys
import csv
from tkinter import Tk, filedialog
import os
import matplotlib.pyplot as plt
import numpy as np

def getbasepeaks(record):

    #FWO_1 = base order only for data 1-4 tho
    #PBAS1 = called base sequence
    #PCON1 = quality score(s)
    abifraw = record.annotations['abif_raw']
    traceA = abifraw['DATA10']
    traceC = abifraw['DATA12']
    traceG = abifraw['DATA9']
    traceT = abifraw['DATA11']
    #ploc 2 is peak indexes where each base call occurs in the traces
    positions = abifraw['PLOC2']
    #actual string dna
    seq = str(record.seq)

    results = []
    #for eahc index of current base, and actual base of sequence at eahc position
    for i, base in enumerate(seq):
        pos = positions[i]
        results.append((base, pos, traceA[pos], traceC[pos], traceG[pos], traceT[pos]))
    return results

def printpeaks(peaks):
    print("Base\tPos\tA\tC\tG\tT")
    baseseq = ""
    for base, pos, a, c, g, t in peaks:
        baseseq += base
        print(f"{base}\t{pos}\t{a}\t{c}\t{g}\t{t}")


def getextended(ab1seq, anchorend):
    newstart = anchorend + 42  #was 41
    newend = anchorend + 79 #45 # was 44
    
    newab1sequence = ab1seq[newstart:newend + 1]
    
    return newstart, newend, newab1sequence

def lookat4bp(ab1seq, startext, endext, seqext, peaks, similaritysequenceCS):
    seqL = len(similaritysequenceCS)
    site = ab1seq[startext:startext + 4]
    print(f" 8th, 9th, 10th, 11th base =  {site}")

    print("\nBase-by-base COmparison @ each index:")
    print("Index (8-11th) = CS Intact Base Expected | Other Bases (NHEJ)  | % CS Intact | % Other | % Calc. NHEJ")
    print("-" * 70)

    results = []
    BaseList =  ['A', 'C', 'G', 'T']

    for i in range(seqL):
        csbase = similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A"
       
        PossibleNHEJBases = []

        for eachbase in BaseList:
            if eachbase != csbase:
                PossibleNHEJBases.append(eachbase)


        if startext + i < len(peaks):
            actualbase, pos, A, C, G, T = peaks[startext + i]
            peakdict = {'A': A, 'C': C, 'G': G, 'T': T}
            totalpeaks = sum(peakdict.values())


            #print(peakdict.get(csbase, 0))
            PercentageCS = (peakdict.get(csbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
        
            OtherCount = 0
            for eachbase2 in PossibleNHEJBases:
                OtherCount += (peakdict.get(eachbase2, 0))
            
            PercentageOther = OtherCount / totalpeaks * 100 if totalpeaks > 0 else 0.0

            NHEJpercentageFINAL = PercentageOther * (4/3) 
  



            



        else:
            #percentageCR = 0.0
            PercentageCS = 0.0
            #difference = 0.0
            PercentageOther = 0.0
            NHEJpercentageFINAL = 0.0


        print(f"{i+1:5} | {csbase:7} | {PossibleNHEJBases} | {PercentageCS:6.1f} | {NHEJpercentageFINAL:8.1f}")
        # ['Index (8-11th)', 'CS Intact Base Expected', 'Other Bases',  '% CS Intact', '% Other' , '% Calc. NHEJ']

        results.append({
            #'CRbase': crbase,
            'CSbase': csbase,
            #'NHEJ1base' : nhej1b,
            #'NHEJ2base' : nhej2b,
            'OtherBases': PossibleNHEJBases,
            #'percentageCR': percentageCR,
            'PercentageCS': PercentageCS,
            'PercentageOther': PercentageOther,
            'PercentageNHEJ': NHEJpercentageFINAL
        })
    


    return results

def main():
    root = Tk()
    root.withdraw()
    
    print("Select AB1 ")
    filenames = filedialog.askopenfilenames(
        title="Select AB1",
        filetypes=[("AB1 files", "*.ab1")]
    )
    if not filenames:
        print("none")
        sys.exit(1)
    
    print("location: ")
    csvpath = filedialog.asksaveasfilename(
        title="Save as",
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv")]
    )
    
    if not csvpath:
        print("not saved.")
        sys.exit(1)

    allresults = []
    samplenames = []
    #similaritysequenceCR = "TATC"
    similaritysequenceCS = "TATCGCCATCCGGGAxTGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGGTGCGCCTATGTCCAG"
    #"CCGG"
    #simseqnhej1 = "AGCT"
    #simseqnhej2 = "GTAA"

    #If not CS, is probably an NHEJ, even if CS, 25% chance it is NHEJ if base is the same. ^^
    
    for filename in filenames:
        print(f"\nProcessing {os.path.basename(filename)}")
        samplenames.append(os.path.basename(filename))
        
        try:
            record = SeqIO.read(filename, "abi")
            peaks = getbasepeaks(record)
            ab1seq = str(record.seq)
            print(ab1seq)
            findpattern = "GTTCCGGTGCCGGAAAGACGACCCT"
            
            anchorstart = ab1seq.find(findpattern)
            if anchorstart == -1:
                raise ValueError(f"Anchor '{findpattern}' not found")
            anchorend = anchorstart + len(findpattern) - 1
    
            startext, endext, seqext = getextended(ab1seq, anchorend)
            percentages = lookat4bp(ab1seq, startext, endext, seqext, peaks, similaritysequenceCS)            
            allresults.append(percentages)
            
        except Exception as e:
            print(f"Error with {filename}: {str(e)}")
            #/N/A 4 bad sample read
            emptyresult = []
            for i in range(len(similaritysequenceCS)):
                emptyresult.append({
                    #'CRbase': "N/A",
                    'CSbase': "N/A",
                    'OtherBases' : "N/A",
                    #'NHEJ1base': "N/A",
                    #'NHEJ2base': "N/A",
                    #'percentageCR': 0.0,
                    'PercentageCS': "N/A",
                    'PercentageOther':"N/A",
                    'PercentageNHEJ': "N/A"
                })
            allresults.append(emptyresult)
            continue
        
    print("\nALL:", samplenames)

    with open(csvpath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for samplename, sampleresult in zip(samplenames, allresults):
            writer.writerow([samplename])
            header = ['Index (8-11th)', 'CS Intact Base Expected', 'Other Bases', '% CS Intact', '% Other' , ' % Calc. NHEJ']
            #print("Index (8-11th) = CS Intact Base Expected | Other Bases  | % CS Intact | % Other | % Calc. NHEJ")

            writer.writerow(header)

            for i in range(len(similaritysequenceCS)):
                otherBstr = ','.join(sampleresult[i]['OtherBases']) if isinstance(sampleresult[i]['OtherBases'], list) else sampleresult[i]['OtherBases']
                row = [i + 1, sampleresult[i]['CSbase'], otherBstr , sampleresult[i]['PercentageCS'], sampleresult[i]['PercentageOther'], sampleresult[i]['PercentageNHEJ']]
                writer.writerow(row)

            writer.writerow([])
    
    print(f"\nSaved @ {csvpath}")
    
    input("done")

if __name__ == "__main__":
    main()