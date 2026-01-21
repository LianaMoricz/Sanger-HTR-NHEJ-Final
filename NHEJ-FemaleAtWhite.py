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
    newstart = anchorend + 39  #was 41
    newend = anchorend + 42 #44
    
    newab1sequence = ab1seq[newstart:newend + 1]
    
    return newstart, newend, newab1sequence

def lookat4bp(ab1seq, startext, endext, seqext, peaks, similaritysequenceCR, similaritysequenceCS, simseqnhej1, simseqnhej2):
    site = ab1seq[startext:startext + 4]
    print(f" 8th, 9th, 10th, 11th base =  {site}")

    print("\nBase-by-base COmparison @ each index:")
    print("Index (8-11th) = CR Base Expected | CS Intact Base Expected | Other (NHEJ)  | % CR | % CS Intact | % Other (NHEJ) | ")
    print("-" * 60)

    results = []
    for i in range(4):
        crbase = similaritysequenceCR[i] if i < len(similaritysequenceCR) else "N/A"
        csbase = similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A"
        nhej1b = simseqnhej1[i] if i < len(simseqnhej1) else "N/A"
        nhej2b = simseqnhej2[i] if i < len(simseqnhej2) else "N/A"




        if startext + i < len(peaks):
            actualbase, pos, A, C, G, T = peaks[startext + i]
            peakdict = {'A': A, 'C': C, 'G': G, 'T': T}
            totalpeaks = sum(peakdict.values())

            percentageCR = (peakdict.get(crbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            percentageCS = (peakdict.get(csbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            percNHEJ1 = (peakdict.get(nhej1b, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            percNHEJ2 = (peakdict.get(nhej2b, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            totalnhej = ((peakdict.get(nhej1b, 0) + (peakdict.get(nhej2b, 0))) / totalpeaks) * 100 if totalpeaks > 0 else 0.0

            NHEJpercentageFINAL = totalnhej * 2

        else:
            percentageCR = 0.0
            percentageCS = 0.0
            difference = 0.0

        print(f"{i+1:5} | {crbase:7} | {csbase:7} | {nhej1b}/{nhej2b:7} | {percentageCR:6.1f} | {percentageCS:6.1f} | {NHEJpercentageFINAL:8.1f}")
        # ['Index (8-11th)',  'CR Base Expected', 'CS Intact Base Expected', 'Other (NHEJ)', '% CR', '% CS Intact', '% Other (NHEJ)']

        results.append({
            'CRbase': crbase,
            'CSbase': csbase,
            'NHEJ1base' : nhej1b,
            'NHEJ2base' : nhej2b,
            'percentageCR': percentageCR,
            'percentageCS': percentageCS,
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
    similaritysequenceCR = "TATC"
    similaritysequenceCS = "CCGG"
    simseqnhej1 = "AGCT"
    simseqnhej2 = "GTAA"
    
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
            percentages = lookat4bp(ab1seq, startext, endext, seqext, peaks, similaritysequenceCR, similaritysequenceCS, simseqnhej1, simseqnhej2)            
            allresults.append(percentages)
            
        except Exception as e:
            print(f"Error with {filename}: {str(e)}")
            #/N/A 4 bad sample read
            emptyresult = []
            for i in range(4):
                emptyresult.append({
                    'CRbase': "N/A",
                    'CSbase': "N/A",
                    'NHEJ1base': "N/A",
                    'NHEJ2base': "N/A",
                    'percentageCR': 0.0,
                    'percentageCS': 0.0,
                    'PercentageNHEJ': 0.0
                })
            allresults.append(emptyresult)
            continue
        
    print("\nALL:", samplenames)

    with open(csvpath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for samplename, sampleresult in zip(samplenames, allresults):
            writer.writerow([samplename])
            header = ['Index (8-11th)', 'CR Base Expected', 'CS Intact Base Expected', 'Other (NHEJ)', '% CR', '% CS Intact', '% Other (NHEJ)']
            writer.writerow(header)

            for i in range(4):
                row = [i + 1, sampleresult[i]['CRbase'], sampleresult[i]['CSbase'], f"{sampleresult[i]['NHEJ1base']}/{sampleresult[i]['NHEJ2base']}", sampleresult[i]['percentageCR'], sampleresult[i]['percentageCS'], sampleresult[i]['PercentageNHEJ']]
                writer.writerow(row)

            writer.writerow([])
    
    print(f"\nSaved @ {csvpath}")
    
    input("done")

if __name__ == "__main__":
    main()