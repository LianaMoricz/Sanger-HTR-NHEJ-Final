from Bio import SeqIO
import sys
import csv
from tkinter import Tk, filedialog
import os

def getbasepeaks(record):
    abifraw = record.annotations['abif_raw']
    traceA = abifraw['DATA10']
    traceC = abifraw['DATA12']
    traceG = abifraw['DATA9']
    traceT = abifraw['DATA11']
    positions = abifraw['PLOC2']
    seq = str(record.seq)

    results = []
    for i, base in enumerate(seq):
        pos = positions[i]
        results.append((base, pos, traceA[pos], traceC[pos], traceG[pos], traceT[pos]))
    return results

def finddeletionregion(ab1seq, anchor="GTTCCGGTGCCGGAAAGACGACCCT"):
    anchorpos = ab1seq.find(anchor)
    if anchorpos == -1:
        raise ValueError("Anchor sequence not found in AB1 file")
    
    deletionstart = anchorpos + len(anchor)
    return deletionstart

def lookat12bpdeletion(ab1seq, startext, endext, peaks, similaritysequenceCR, similaritysequenceCS):
    firstdeletionsite = ab1seq[startext:startext + 12]
    print(f"Deletion Site =  {firstdeletionsite}")

    print("\nBase-by-base CR vs CS :")
    print("Index | CR Base | CS Base | %CR | %CS | Difference")
    print("-" * 60)

    results = []
    for i in range(12):
        crbase = similaritysequenceCR[i] if i < len(similaritysequenceCR) else "N/A"
        csbase = similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A"

        if startext + i < len(peaks):
            actualbase, pos, A, C, G, T = peaks[startext + i]
            peakdict = {'A': A, 'C': C, 'G': G, 'T': T}
            totalpeaks = sum(peakdict.values())

            percentageCR = (peakdict.get(crbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            percentageCS = (peakdict.get(csbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0

            if percentageCR + percentageCS > 0:
                #difference = ((percentageCR - percentageCS) / (percentageCR + percentageCS)) * 100
                # = with ATG ((percentageCR) / (percentageCR + percentageCS)) * 100
                difference = ((percentageCR) / (percentageCR + percentageCS)) * 100
            else:
                difference = 0.0

        else:
            percentageCR = 0.0
            percentageCS = 0.0
            difference = 0.0

        print(f"{i+1:5} | {crbase:7} | {csbase:7} | {percentageCR:6.1f} | {percentageCS:6.1f} | {difference:8.1f}")
        results.append({
            'CRbase': crbase,
            'CSbase': csbase,
            'percentageCR': percentageCR,
            'percentageCS': percentageCS,
            'difference': difference
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
    
    print("loaction: ")
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
    similaritysequenceCR = "GCTGAATGCCCT"
    similaritysequenceCS = "TGCCTTTCGATC"
    
    for filename in filenames:
        print(f"\nProcessing {os.path.basename(filename)}")
        samplenames.append(os.path.basename(filename))
        
        try:
            record = SeqIO.read(filename, "abi")
            peaks = getbasepeaks(record)
            ab1seq = str(record.seq)
            
            #Findwondowstart
            windowstart = finddeletionregion(ab1seq)
            print(f"Found deletion region starting at position: {windowstart}")
            
            #12bp deletion region
            percentages = lookat12bpdeletion(
                ab1seq, 
                windowstart, 
                windowstart + 11,  # 12
                peaks, 
                similaritysequenceCR, 
                similaritysequenceCS
            )            
            allresults.append(percentages)
            
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            allresults.append([{'CRbase': similaritysequenceCR[i] if i < len(similaritysequenceCR) else "N/A", 
                              'CSbase': similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A", 
                              'percentageCR': 0.0, 'percentageCS': 0.0, 'difference': 0.0} for i in range(12)])
            continue
        
    print("\nAll samples = ", samplenames)

    with open(csvpath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for samplename, sampleresult in zip(samplenames, allresults):
            writer.writerow([samplename])
            header = ['Index', 'CR Base', 'CS Base', 'Percentage CR', 'Percentage CS', 'Difference']
            writer.writerow(header)

            for i in range(12):
                row = [i + 1, sampleresult[i]['CRbase'], sampleresult[i]['CSbase'],
                    sampleresult[i]['percentageCR'], sampleresult[i]['percentageCS'], sampleresult[i]['difference']]
                writer.writerow(row)

            writer.writerow([])
    
    print(f"\nSaved @ {csvpath}")
    
    input("done")

if __name__ == "__main__":
    main()