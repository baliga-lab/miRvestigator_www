import sys, re, os, math, shutil
from subprocess import *
from copy import deepcopy
from random import sample, randint
import cPickle

# Custom libraries
from miRvestigator import miRvestigator
from pssm import pssm

# Library for Pyro - Remote Object Communication
import Pyro.core

# Reverse complement
def reverseMe(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(seq)

# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
# The weeder program can be found at:  http://159.149.109.9/modtools/
# I modified the C code and recompiled to make Weeder look for the FreqFiles
# folder in /local/FreqFiles. Then I made symbolic links in my PATH so that
# weeder could be run from the command line as weederlauncher. You will also
# have to add weederTFBS.out and adviser.out to the PATH in order to run.
def weeder(seqFile=None, percTargets=50, revComp=False, bgModel='HS'):
    if not os.path.exists('tmp/weeder'):
        os.makedirs('tmp/weeder')
    
    # First run weederTFBS for 6bp motifs
    weederArgs = ' '+str(seqFile)+' '+str(bgModel)+' small T50'
    if revComp==True:
        weederArgs += ' -S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weederlauncher " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    """# First run weederTFBS for 6bp motifs
    weederArgs = '-f '+str(seqFile)+' -W 6 -e 1 -O HS -R '+str(percTargets)
    if revComp==True:
        weederArgs += ' -S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weeder " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    # Second run weederTFBS for 8bp motifs
    weederArgs = '-f '+str(seqFile)+' -W 8 -e 2 -O HS -R '+str(percTargets)
    if revComp==True:
        weederArgs += ' -S'
    weederProc = Popen("weeder " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    
    # Finally run adviser
    weederArgs = str(seqFile)
    weederProc = Popen("adviser " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()
    errOut.close()
    """
    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqFile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break
    
    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    while 1:
        name = outLines.pop(0).strip() # Get match
        if not name.find('(not highest-ranking)') == -1:
            break
        # Get redundant motifs
        outLines.pop(0)
        redMotifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [i for i in col.strip().split('\t')[1].split(' ') if i]
            colSum = 0
            for i in nums:
                colSum += int(i.strip())
            matrix += [[ float(nums[0])/float(colSum), float(nums[1])/float(colSum), float(nums[2])/float(colSum), float(nums[3])/float(colSum)]]
            col = outLines.pop(0)
        PSSMs += [pssm(biclusterName=name,nsites=instances,eValue=hitBp[len(matrix)][1],pssm=matrix,genes=redMotifs)]
    return PSSMs


# Reverse complement
def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

# Complement
def complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'U':'A'}
    complseq = [complement[base] for base in seq]
    return complseq

# Convert to RNA
def conv2rna(seq):
    conversion = {'A':'A', 'T':'U', 'C':'C', 'G':'G', 'N':'N', 'U':'U'}
    rnaSeq = [conversion[base] for base in list(seq)]
    return ''.join(rnaSeq)


def alignSeed(alignment, seed, motif):
    alignment.pop() # Get rid of the extra state which is added by the forwardViterbi function
    start = 1    
    if alignment[0]=='NM1':
        for i in alignment:
            if i=='NM1':
                start += 1
    # Alignment
    motifAlign = '<font color=\'#CC0000\'><b>Motif</b></font><font color=\'#ffffff\'>_</font>5\'<font color=\'#ffffff\'>_</font>'
    seed = list(conv2rna(complement(seed)))
    seedAlign = '<font color=\'#ffffff\'>______</font>3\'<font color=\'#ffffff\'>_</font>'
    motif = list(conv2rna(motif))
    alignMe = alignment
    aligned = '<font color=\'#ffffff\'>_________</font>'
    lenMatch = 0
    # First get them zero'd to the same point
    if start>1:
        for i in range(start-1):
            seedAlign += seed.pop(0)
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += '-'
            alignMe.pop(0)
    if len(alignMe)>0 and not alignMe[0]=='PSSM0' and not alignMe[0]=='WOBBLE0':
        if alignMe[0][0]=='P':
            upTo = int(alignMe[0][4])
        elif alignMe[0][0]=='W':
            upTo = int(alignMe[0][6])
        for i in range(upTo):
            seedAlign += '<font color=\'#cccccc\'>-</font>'
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += motif.pop(0)
    # Then add the parts that align
    while 1:
        if len(alignMe)==0 or alignMe[0]=='NM2':
            break
        seedAlign += seed.pop(0)
        if alignMe[0][0]=='P':
            aligned += '<font color=\'#ff0000\'>|</font>'
        elif alignMe[0][0]=='W':
            aligned += '<font color=\'#0000ff\'>:</font>'
        lenMatch += 1
        motifAlign += motif.pop(0)
        alignMe.pop(0)
    # Then do the ending
    if len(alignMe)>0:
        for i in alignMe:
            seedAlign += seed[0]
            seed = seed[1:]
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += '<font color=\'#cccccc\'>-</font>'
        alignMe = []
    if len(motif)>0 and len(alignMe)==0:
        for i in motif:
            seedAlign += '<font color=\'#cccccc\'>-</font>'
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += i
    motifAlign += '<font color=\'#ffffff\'>_</font>3\'<font color=\'#ffffff\'>___________</font>'
    aligned +=    '<font color=\'#ffffff\'>______________</font>'
    seedAlign +=  '<font color=\'#ffffff\'>_</font>5\'<font color=\'#ffffff\'>_</font><font color=\'#cc0000\'><b>miRNA Seed</b></font>'
    return [motifAlign, aligned, seedAlign, lenMatch]

class miRwww(Pyro.core.ObjBase):
    def __init__(self):
        Pyro.core.ObjBase.__init__(self)
    def run(self, genes, seedModels, wobble, cut, bgModel, motifSizes):
        cut = float(cut)
        curRunNum = randint(0,1000000)

        # 1. Read in sequences
        seqFile = open('p3utrSeqs_Homo_sapiens.csv','r')
        seqLines = seqFile.readlines()
        ids = [i.strip().split(',')[0].upper() for i in seqLines]
        sequences = [i.strip().split(',')[1] for i in seqLines]
        seqs = dict(zip(ids,sequences))
        seqFile.close()

        # 2. Get sequences for each target
        miRSeqs = {}
        for gene in genes:
            if gene in seqs:
                miRSeqs[gene] = seqs[gene]

        # 3. Make a FASTA file
        if not os.path.exists('tmp/fasta'):
            os.makedirs('tmp/fasta')
        fastaFile = open('tmp/fasta/tmp'+str(curRunNum)+'.fasta','w')
        for seq in miRSeqs:
            fastaFile.write('>'+str(seq)+'\n'+str(miRSeqs[seq])+'\n')
        fastaFile.close()
        
        # 4. Run weeder
        weederPSSMs1 = weeder(seqFile='tmp/fasta/tmp'+str(curRunNum)+'.fasta', percTargets=50, revComp=False, bgModel=bgModel)
        
        # 4a. Take only selected size motifs
        weederPSSMsTmp = []
        for pssm1 in weederPSSMs1:
            if 6 in motifSizes and len(pssm1.getName())==6:
                weederPSSMsTmp.append(deepcopy(pssm1))
            if 8 in motifSizes and len(pssm1.getName())==8:
                weederPSSMsTmp.append(deepcopy(pssm1))
        weederPSSMs1 = deepcopy(weederPSSMsTmp)
        del weederPSSMsTmp

        # 5. Run miRvestigator HMM
        mV = miRvestigator(weederPSSMs1,seqs.values(),seedModel=seedModels,minor=True,p5=True,p3=True,wobble=wobble,wobbleCut=cut,textOut=False)
        
        # 6. Read in miRNAs to get mature miRNA ids
        import gzip
        miRNAFile = gzip.open('mature.fa.gz','r')
        miRNADict = {}
        while 1:
            miRNALine = miRNAFile.readline()
            seqLine = miRNAFile.readline()
            if not miRNALine:
                break
            # Get the miRNA name
            miRNAData = miRNALine.lstrip('>').split(' ')
            curMiRNA = miRNAData[0]
            if (curMiRNA.split('-'))[0]=='hsa':
                miRNADict[curMiRNA] = miRNAData[1]
        miRNAFile.close()

        # 6. Clean-up after yerself
        os.remove('tmp/fasta/tmp'+str(curRunNum)+'.fasta')
        os.remove('tmp/fasta/tmp'+str(curRunNum)+'.fasta.wee')
        os.remove('tmp/fasta/tmp'+str(curRunNum)+'.fasta.mix')
        os.remove('tmp/fasta/tmp'+str(curRunNum)+'.fasta.html')

        # 7. Return results
        s = '<html><body bgcolor=\'#333333\' link=\'cc0000\' vlink=\'cc0000\'><font face=\'arial\'><center><table bgcolor=\'#999966\' cellpadding=\'10%\'><tr><td><center>'
        for pssm1 in weederPSSMs1:
            s += '<table width=\'100%\' bgcolor=\'#333333\' cellpadding=\'15%\'><tr><td align=\'center\' valign=\'center\'><font size=6><b><font color=\'#ff0000\'>miR</font><font color=\'#cccc00\'>vestigator Framework Results</font></b></font></td></tr></table>'
            s += '<p><table width=\'100%\' bgcolor=\'#333333\' cellpadding=\'15%\'><tr><td align=\'center\' valign=\'center\'><font size=4><b><font color=\'#cccc00\'>Top Ten miRNAS Matching the Weeder Motif</font> <font color=\'#ff0000\'>'+str(pssm1.getName())+'</font></b></font></td></tr></table>'
            scoreList = mV.getScoreList(pssm1.getName())
            s += '<table width=\'100%\' cellpadding=\'15%\'><tr><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>miRNA Name</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>miRNA Seed</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Seed Model</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Length of</br>Alignment</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Alignment</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Viterbi P-Value</font></b></center></td></tr>'
            for k in range(10):
                i = scoreList[k]
                align1 = alignSeed(i['statePath'], i['miRNA.seed'], pssm1.getName())
                s += '<tr><td bgcolor=\'#ffffff\'><center>'+str('</br>'.join(['<a href=\'http://mirbase.org/cgi-bin/mature.pl?mature_acc='+str(miRNADict[j.strip()])+'\' target=\'_blank\'>'+str(j.strip())+'</a>' for j in i['miRNA.name'].split('_')]))+'</center></td><td bgcolor=\'#ffffff\'><center>'+conv2rna(reverseComplement(str(i['miRNA.seed'])))+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(i['model'])+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(align1[3])+'</center></td><td bgcolor=\'#ffffff\'>'
                s += '<center><pre>'+str(align1[0])+'\n'+str(align1[1])+'\n'+str(align1[2])+'</pre></center>'
                s += '</td><td bgcolor=\'#ffffff\'><center>'+str('%.2g' % float(i['vitPValue']))+'</center></td></tr>'
            s += '</table></p>'
            #'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')')
            # pssm1.nsites
            s += '<p><table width=\'100%\' bgcolor=\'#333333\' cellpadding=\'15%\'><tr><td align=\'center\' valign=\'center\'><font size=4><b><font color=\'#cccc00\'>Position of Putative miRNA Binding Sites in Submitted Genes</br>for the Weeder Motif</font> <font color=\'#ff0000\'>'+str(pssm1.getName())+'</font></b></font></td></tr></table>'
            s += '<table width=\'100%\' cellpadding=\'15%\'><tr><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Entrez Gene ID</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Site</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Start Relative to</br>Stop Codon</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>% Similarity to Consensus Motif</br>(</font><font color=\'#cc0000\'>High Quality</font><font color=\'#ffffff\'> | </font><font color=\'#cccc00\'>Medium Quality</font><font color=\'#ffffff\'> | </font><font color=\'#00cc00\'>Fair Quality</font><font color=\'#ffffff\'>)</font></b></center></td></tr>'
            for i in pssm1.nsites:
                col1 = '#000000'
                if float(i['match']) >= float(95):
                    col1 = '#cc0000'
                elif float(i['match']) >= float(90):
                    col1 = '#cccc00'
                elif float(i['match']) >= float(85):
                    col1 = '#00cc00'
                s += '<tr><td bgcolor=\'#ffffff\'><center>'+str('<a href=\'http://www.ncbi.nlm.nih.gov/gene/'+str(i['gene'])+'\' target=\'_blank\'>'+str(i['gene'])+'</a>')+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(i['site'])+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(i['start'])+'</center></td><td bgcolor=\'#ffffff\'><font color=\''+str(col1)+'\'><center><b>'+i['match']+'</b></center></font></td></tr>'
            s += '</table></p>'
        s += '</center></td></tr></table></center></font></body></html>'
        return s

Pyro.core.initServer()
daemon = Pyro.core.Daemon()
miRwww = miRwww()
uri = daemon.connect(miRwww, 'miRwww')

print uri
uriOut = open('/var/www/uri','w')
uriOut.write(str(uri))
uriOut.close()

daemon.requestLoop()

