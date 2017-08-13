__author__ = 'George'


import math
import os,sys
from multiprocessing import Process,cpu_count
from multiprocessing import Pool

try:  import psutil
except ImportError:
    # print('Warning! psutil package is recommended but not found')
    pass

try:  import resource
except ImportError:
    # print('Warning! \'resource\' package is recommended but not found')
    pass



def CorrectThreadCountAccordingToAvailableMem(ThreadsCount):
    available_mem_gib = None
    try:
        available_mem_gib = psutil.phymem_usage().total / float (2 ** 30) * (100 - psutil.phymem_usage().percent) / 100
    except:
        pass

    used_mem_gib = None
    try:
        process = psutil.Process(os.getpid())
        used_mem_gib = process.get_memory_info()[0] / float(2 ** 30)
    except:
        try:
            used_mem_gib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / float(2**20)
        except: pass


    if used_mem_gib != None and available_mem_gib != None:
        print('memory used: %g Gb'%used_mem_gib)
        print('memory available: %g Gb'%available_mem_gib)
        ThreadsCount = min(ThreadsCount,4+int(available_mem_gib/used_mem_gib*5))
        print('Thread count was corrected to %d'%ThreadsCount)
        return ThreadsCount
    else:
        print('Cannot identify sizes of total, available and used memory. Threads count will not be corrected')
        return ThreadsCount


def Omniconvert(x):
    try: return int(x)
    except:
        try: return float(x)
        except: return x

def ReadTableFromTXT_Autoformat(FileName,delimiter='\t',skiprows=0,comments='#'):
    print('Loading %s...'%FileName)
    F = open(FileName,'r',buffering = 100000000)
    for row in range(skiprows): F.readline()
    FileSize = os.path.getsize(FileName)
    ReadenSize = 0
    Chunk = 0
    while True:
        L = F.readline()
        ReadenSize += len(L)
        if Chunk > 2000:
            sys.stdout.write('\rCompleted %.2f%%'%(100*ReadenSize/FileSize))
            Chunk = 0
        Chunk += 1
        if len(L) == 0: break
        if L.startswith(comments): continue
        yield [Omniconvert(y) for y in L.replace('\r','').replace('\n','').split(delimiter)]
    print('\tCompleted.                       ')
    F.close()

def LoadTableFromTXT_Autoformat(FileName,delimiter='\t',skiprows=0,comments='#'):
    F = open(FileName,'r',buffering = 100000000)
    if comments != None:
        Table = [[Omniconvert(y) for y in x.replace('\r','').replace('\n','').split(delimiter)] for x in F.readlines()[skiprows:] if not x.startswith(comments)]
    else:
        Table = [[Omniconvert(y) for y in x.replace('\r','').replace('\n','').split(delimiter)] for x in F.readlines()[skiprows:]]
    F.close()
    return Table


class T_ENCODE_Gene():
    def __init__(self,ID,Chr,Strand):
        self.ID = ID
        self.Chr = Chr
        self.Strand = Strand
        self.Name = None
        self.Transcripts = []
        self.Proteins = []
        self.Start = None
        self.End = None

        self.Jaspar_Hits_by_TF_name = {}

        self.Current_LookupDistances_upstream = []
        self.Current_LookupDistances_downstream = []
        self.TSSs = []
        self.ENCODE_TFBS_by_TSS_Number = []
        self.SegDescriptions_by_TSS_Number = []

    def CalcStartEnd(self):
        self.Start = 10000000000
        self.End = 0
        for T in self.Transcripts:
            for (EStart,EEnd) in zip(T.ExonStarts,T.ExonEnds):
                self.Start = min(EStart,self.Start)
                self.End = max(EEnd,self.End)

class T_ENCODE_Transcript():
    def __init__(self,ID,Chr,Strand,GeneID):
        self.ID = ID
        self.Name = None
        self.Chr = Chr
        self.GeneID = GeneID
        self.GeneName = None
        self.ProteinID = None
        self.ExonStarts = []
        self.ExonEnds = []
        self.Strand = Strand

class T_ENCODE_Protein():
    def __init__(self,ID,Chr,Strand,GeneID,TranscriptID):
        self.ID = ID
        self.Name = None
        self.GeneID = GeneID
        self.TranscriptID = TranscriptID
        self.Chr = Chr
        self.ExonStarts = []
        self.ExonEnds = []
        self.Strand= Strand

def gmean(Numbers):
    x = 1
    for N in Numbers:  x *= N
    return x**(1/len(Numbers))

def Load_GTF_for_ENCODE(FileName, ChrLimit=None):
    TranscriptsByID = {}
    ProteinsByID = {}
    GenesByID = {}
    GenesByName = {}

    InGtfFile = open(FileName,'r',buffering=10000000)
    print('Loading GTF file...')

    SourceFileSize=os.path.getsize(FileName)
    LineCount = 0
    Chunk = 0
    Line=InGtfFile.readline()
    while Line!='':
        LineCount+=1
        Chunk=Chunk+1
        if Chunk>9999:
            percents=InGtfFile.tell()/SourceFileSize*100
            sys.stdout.write("\rCompleted:%.1f%%. Processed %d strings" %(percents,LineCount))
            Chunk=0

        if Line[0]=='#':
            Line=InGtfFile.readline()
            continue
        if Line[-1].isprintable()==False:
            Line=Line[:-1]
        if Line[-1].isprintable()==False:
            Line=Line[:-1]

        Cells=Line.split('\t')
        Chr=Cells[0]
        if ChrLimit!=None:
            if not Chr in ChrLimit:
                Line=InGtfFile.readline()
                continue

        if Cells[1] in ['nonsense_mediated_decay','antisense','retained_intron','processed_transcript']:
            Line=InGtfFile.readline()
            continue


        In=int(Cells[3])-1
        Out=int(Cells[4])-1
        Strand = '-' if (Cells[6]=='-') else '+'
        Des=(Cells[-1].replace(' ','').replace(';;',';').split(';'))[:-1]
        Items={}
        for D in Des:
            Items[D.split('"')[0]]=D.split('"')[1]

        if Cells[2]=='exon':
            try:
                transcript_id = Items['transcript_id']
                gene_id = Items['gene_id']
            except KeyError:
                print('Incorrect GTF file: no transcript_id or gene_id field for EXON. Line %s'%Line)
                exit()

            if not transcript_id in TranscriptsByID.keys():
                TranscriptsByID[transcript_id] = T_ENCODE_Transcript(transcript_id,Chr,Strand,gene_id)

            TranscriptsByID[transcript_id].ExonStarts.append(In)
            TranscriptsByID[transcript_id].ExonEnds.append(Out)

            if 'gene_name' in Items: TranscriptsByID[transcript_id].GeneName = Items['gene_name']

        elif Cells[2]=='CDS':
            try:
                protein_id = Items['protein_id']
                transcript_id = Items['transcript_id']
                gene_id = Items['gene_id']
            except KeyError:
                print('Incorrect GTF file: no transcript_id, protein_id or gene_id field for CDS. Line %s'%Line)
                exit()

            if not protein_id in ProteinsByID.keys():
                ProteinsByID[protein_id] = T_ENCODE_Protein(protein_id,Chr,Strand,gene_id,transcript_id)

            ProteinsByID[protein_id].ExonStarts.append(In)
            ProteinsByID[protein_id].ExonEnds.append(Out)
        Line=InGtfFile.readline()

    print('\rCompleted. Distributing genes, transcripts and proteins....')
    for transcript_id in sorted(TranscriptsByID.keys()):
        T = TranscriptsByID[transcript_id]
        if not T.GeneID in GenesByID.keys():
            GenesByID[T.GeneID] = T_ENCODE_Gene(T.GeneID,T.Chr,T.Strand)
            if T.GeneName != None:
                GenesByID[T.GeneID].Name = T.GeneName
        GenesByID[T.GeneID].Transcripts.append(T)
    for protein_id in sorted(ProteinsByID.keys()):
        P = ProteinsByID[protein_id]
        if not P.GeneID in GenesByID.keys():
            GenesByID[P.GeneID] = T_ENCODE_Gene(P.GeneID,P.Chr,P.Strand)
            print('Warning: protein %s with gene_id %s is present whereas corresponding transcript with the same gene_id lacks'%(protein_id,P.GeneID))
        GenesByID[P.GeneID].Proteins.append(P)

        if not P.TranscriptID in TranscriptsByID.keys():
            print('Warning: protein %s with transcript_id %s is present whereas this transcript lacks'%(protein_id,P.TranscriptID))

    for G in GenesByID.values():
        for T in G.Transcripts:
            T.ExonStarts.sort()
            T.ExonEnds.sort()
        if G.Name != None:
            GenesByName[G.Name] = G
        G.CalcStartEnd()

    print('Completed. Found %d genes, %d transcripts, %d proteins'%(len(GenesByID),len(TranscriptsByID),len(ProteinsByID)))

    return GenesByName,GenesByID,TranscriptsByID,ProteinsByID

def Load_Genome(GenomeFileName,ChrLimit = None):
    print('Loading genome file....\n')
    ChrSeqs = {}
    GenomeFile=open(GenomeFileName,'r',buffering=100000000)
    for record in SeqIO.parse(GenomeFile, "fasta"):
        if ChrLimit!=None:
            if not record.id in ChrLimit:
                continue
        ChrSeqs[record.id]=record
        sys.stdout.write("\rProcessed record '%s' with length %.1f mb    " %(record.id,len(record)/1e+6))

    sys.stdout.write("\rCompleted: 100%%. Processed %d entries              "%len(ChrSeqs))
    GenomeFile.close()
    return ChrSeqs


class TGenomeSegmentationDescription():
    def __init__(self,FileName,CellType,ML_Type):
        self.CellType = CellType
        self.ML_Type = ML_Type
        self.FileName = FileName
        self.SegListsByChr = {}

class TGenomeSegment():
    def __init__(self,Start,End,Type):
        self.Start = Start
        self.End = End
        self.Type = Type

def LoadGenomeSegmentsFromFile(FileName):
    GenomeSegmentsByChr = {}
    Table = LoadTableFromTXT_Autoformat(FileName)
    for El in Table:
        Chr = El[1][3:] if El[1].casefold().startswith('chr') else El[1]
        Start = El[2]
        End = El[3]
        SegType = El[4]
        if not (type(Start) is int and type(End) is int):
            print('Incorrect Genome segmentation file ' + FileName)
            print('string ' + '\t'.join(El))
            exit()
        Seg = TGenomeSegment(Start,End,SegType)
        try: GenomeSegmentsByChr[Chr].append(Seg)
        except KeyError: GenomeSegmentsByChr[Chr] = [Seg]

    for Chr in sorted(GenomeSegmentsByChr.keys()):
        GenomeSegmentsByChr[Chr].sort(key=(lambda x: x.Start))

    return GenomeSegmentsByChr

def RevealSegTypeByChrByPosition(GenomeSegmentsByChr,Chr,Pos):
    if not Chr in GenomeSegmentsByChr.keys():
        return None
    List = GenomeSegmentsByChr[Chr]
    StartArrayPos = 0
    EndArrayPos = len(List)-1
    while True:
        if List[StartArrayPos].Start <= Pos and List[StartArrayPos].End >= Pos: return List[StartArrayPos].Type
        if List[EndArrayPos].Start <= Pos and List[EndArrayPos].End >= Pos: return List[EndArrayPos].Type
        if EndArrayPos - StartArrayPos == 1: return None

        if List[StartArrayPos].Start > Pos: return None
        if List[EndArrayPos].Start <= Pos:
            if List[EndArrayPos].End >= Pos: return List[EndArrayPos].Type
            else: return None

        MidArrayPos = int((StartArrayPos + EndArrayPos)/2)
        if List[MidArrayPos].Start <= Pos and List[MidArrayPos].End >= Pos: return List[MidArrayPos].Type
        if List[MidArrayPos].Start <= Pos:
            StartArrayPos = MidArrayPos
            continue
        else:
            EndArrayPos = MidArrayPos

def RevealSegIndexInListByPosition(GenomeSegments,Pos,TendTo='end'):
    List = GenomeSegments
    StartArrayPos = 0
    EndArrayPos = len(List)-1
    while True:
        if List[StartArrayPos].Start <= Pos and List[StartArrayPos].End >= Pos: return StartArrayPos
        if List[EndArrayPos].Start <= Pos and List[EndArrayPos].End >= Pos: return EndArrayPos
        if EndArrayPos - StartArrayPos == 1:
            if TendTo == 'end': return EndArrayPos
            return StartArrayPos

        if List[StartArrayPos].Start > Pos: return None
        if List[EndArrayPos].Start <= Pos:
            if List[EndArrayPos].End >= Pos: return EndArrayPos
            else: return None

        MidArrayPos = int((StartArrayPos + EndArrayPos)/2)
        if List[MidArrayPos].Start <= Pos and List[MidArrayPos].End >= Pos: return MidArrayPos
        if List[MidArrayPos].Start <= Pos:
            StartArrayPos = MidArrayPos
            continue
        else:
            EndArrayPos = MidArrayPos


def GetSegDescriptionsForSpecificArea(SegDescriptions,Chr,Start,End,refine = False):
    SegDescriptions_for_specific_area = []
    for SegDes in SegDescriptions:
        SegDes_for_specific_area = TGenomeSegmentationDescription(SegDes.FileName,SegDes.CellType,SegDes.ML_Type)
        SegDescriptions_for_specific_area.append(SegDes_for_specific_area)
        if not Chr in SegDes.SegListsByChr.keys(): continue
        StartIndex = RevealSegIndexInListByPosition(SegDes.SegListsByChr[Chr],Start,TendTo='start')
        if StartIndex is None: StartIndex = 0
        else: StartIndex = max(0,StartIndex-2)

        EndIndex = RevealSegIndexInListByPosition(SegDes.SegListsByChr[Chr],End,TendTo='end')
        if EndIndex is None: EndIndex = max(len(SegDes.SegListsByChr[Chr]) - 1, 0)
        else: EndIndex = min(len(SegDes.SegListsByChr[Chr]) - 1, EndIndex + 2)
        if refine:
            SegDes_for_specific_area.SegListsByChr[Chr] = []
            for index in range(StartIndex,EndIndex+1):
                if SegDes.SegListsByChr[Chr][index].Start <= End and SegDes.SegListsByChr[Chr][index].End >= Start:
                    SegDes_for_specific_area.SegListsByChr[Chr].append(SegDes.SegListsByChr[Chr][index])
        else:
            SegDes_for_specific_area.SegListsByChr[Chr] = SegDes.SegListsByChr[Chr][StartIndex : EndIndex+1]

    return SegDescriptions_for_specific_area


class TENCODE_ChIPSeq_TFBS():
    def __init__(self,TF_Name,Chr,Start,End,Score,ExperimentsCount):
        self.TF_Name = TF_Name
        self.Chr = Chr
        self.Start = Start
        self.End = End
        self.Score = Score
        self.ExperimentsCount = ExperimentsCount

def Get_ENCODE_ChIPSeq_TFBS_for_specified_interval(ENCODE_ChIPSeq_TFBS_by_Chr,Chr,Start,End,Flashback_distance = 50000):
    # Flashback_distance is a distance (bp) in which TFBS are looked upstream in
    # ENCODE_ChIPSeq_TFBS_by_Chr[Chr] list when Start position is found in specificed interval
    # !!!!ENSURE THAT ENCODE_ChIPSeq_TFBS_by_Chr[Chr] is SORTED SOR EACH CHROMOSOME!!!!!
    if not Chr in ENCODE_ChIPSeq_TFBS_by_Chr.keys(): return []
    TFSB_list = ENCODE_ChIPSeq_TFBS_by_Chr[Chr]
    if TFSB_list == []: return {}
    From = 0
    To = len(TFSB_list) - 1
    while True:
        if To - From <= 1: break
        if TFSB_list[To].Start == TFSB_list[From].Start: break
        if TFSB_list[From].Start >= Start:
            To = From
            break
        if TFSB_list[To].Start <= Start:
            From = To
            break
        Mid = int((From + To)/2)
        if TFSB_list[Mid].Start >= Start:
            To = Mid
        else:
            From = Mid
    
    Indexes_of_dropped_items_in_primary_TFSB_list = []
    for index in range(From,To + 1):
        if index > 0 and index < len(TFSB_list):
            if TFSB_list[index].Start <= End and TFSB_list[index].End >= Start:
                Indexes_of_dropped_items_in_primary_TFSB_list.append(index)
    
    ### flashback
    index = From - 1
    while True:
        if index > 0 and index < len(TFSB_list):
            if TFSB_list[index].Start <= End and TFSB_list[index].End >= Start:
                Indexes_of_dropped_items_in_primary_TFSB_list.append(index)
            if Start - TFSB_list[index].Start > Flashback_distance: break
        elif index <= 0: break
        index -= 1

    ### "flashforward" to guarantee best results
    index = To + 1
    while True:
        if index > 0 and index < len(TFSB_list):
            if TFSB_list[index].Start <= End and TFSB_list[index].End >= Start:
                Indexes_of_dropped_items_in_primary_TFSB_list.append(index)
            if End < TFSB_list[index].Start: break
        elif index >= len(TFSB_list) - 1: break
        index += 1
    
    Indexes_of_dropped_items_in_primary_TFSB_list = sorted(set(Indexes_of_dropped_items_in_primary_TFSB_list))
    return [TFSB_list[x] for x in Indexes_of_dropped_items_in_primary_TFSB_list]

        
def AdjustScoreAccordingToTSS_dist_ENCODE(Score,TSS_dist):
    if TSS_dist > 0:
        if TSS_dist < 300:  multiplier = 2.8
        elif TSS_dist < 500:  multiplier = 2.5
        elif TSS_dist < 700: multiplier = 2.4
        elif TSS_dist < 1200: multiplier = 1.7
        elif TSS_dist < 2000: multiplier = 1.2
        elif TSS_dist < 3000: multiplier = 0.8
        else: multiplier = 0.5
    else:
        if TSS_dist > -500:  multiplier = 2.8
        elif TSS_dist > -1000:  multiplier = 2.5
        elif TSS_dist > -2000: multiplier = 2.3
        elif TSS_dist > -5000: multiplier = 1.8
        elif TSS_dist > -10000: multiplier = 1.3
        elif TSS_dist > -20000: multiplier = 0.8
        else: multiplier = 0.5
    return multiplier*Score

class TENCODE_ChIPSeq_Gene_TF_Features():
    def __init__(self):
        self.AdjustedScore = 0
        self.Score = 0
        self.ExperimentsCount = 0
        self.Genome_States_Count_by_Type = {}

def SupplyGenesWith_segment_and_TFBS_info_from_ENCODE(GenesByName,ENCODE_ChIPSeq_TFBS_by_Chr,Flashback_distance,LookupDistance_upstream,LookupDistance_downstream,SegDescriptions=[]):
    print('Supplying genes with ENCODE info...')
    GeneN = 0
    for GeneName in sorted(GenesByName.keys()):
        G = GenesByName[GeneName]
        sys.stdout.write('\rProc.gene %s, %d of %d...'%(G.Name,GeneN+1,len(GenesByName)))
        G.TSSs = []
        if G.Strand == '+':
            G.TSSs = sorted(set([T.ExonStarts[0] for T in G.Transcripts]))
        elif G.Strand == '-':
            G.TSSs = sorted(set([T.ExonEnds[-1] for T in G.Transcripts]))
        else:
            print('[8] Smth strange. Gene %s. Maybe incorrect GTF file.'%(G.Name))
        G.Current_LookupDistances_upstream = []
        G.Current_LookupDistances_downstream = []
        G.ENCODE_TFBS_by_TSS_Number = []
        G.SegDescriptions_by_TSS_Number = []
        for TSS in G.TSSs:
            if G.Strand == '+':
                G.Current_LookupDistances_upstream.append(min(LookupDistance_upstream,TSS))
                G.Current_LookupDistances_downstream.append(LookupDistance_downstream)
            else:
                G.Current_LookupDistances_upstream.append(min(LookupDistance_downstream,TSS))
                G.Current_LookupDistances_downstream.append(LookupDistance_upstream)
            MinCoord = TSS - G.Current_LookupDistances_upstream[-1]
            MaxCoord = TSS + G.Current_LookupDistances_downstream[-1]
            G.ENCODE_TFBS_by_TSS_Number.append(Get_ENCODE_ChIPSeq_TFBS_for_specified_interval(ENCODE_ChIPSeq_TFBS_by_Chr,G.Chr,MinCoord,MaxCoord,Flashback_distance))

            if SegDescriptions != []:
                G.SegDescriptions_by_TSS_Number.append(GetSegDescriptionsForSpecificArea(SegDescriptions,G.Chr,MinCoord,MaxCoord,refine = False))
        GeneN += 1
    print('\rCompleted.                          ')


class Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name_Multiprocess_StartupPars():
    def __init__(self,ThreadN,GenesByName):
        self.GenesByName = GenesByName
        self.ThreadN = ThreadN

def Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name_Multiprocess(Pars):
    print('\r[Thread %d started]. Performing ENCODE ChIP-Seq TFBS assignment to genes and calculating scores...'%ThreadN)
    GenesByName = Pars.GenesByName
    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name = {}
    GeneN = 0
    for GeneName in GenesByName.keys():
        G = GenesByName[GeneName]
        sys.stdout.write('\r[Thread %d] Proc.gene %s, %d of %d...'%(Pars.ThreadN,G.Name,GeneN+1,len(GenesByName)))

        # MinCoord = None
        # MaxCoord = None
        Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position = {}
        for TSS_Number in range(len(G.TSSs)):
            TSS = G.TSSs[TSS_Number]
            # Current_LookupDistance_upstream = Current_LookupDistances_upstream[TSS_Number]
            # Current_LookupDistance_downstream =  Current_LookupDistances_downstream[TSS_Number]
            Current_ENCODE_TFBS = G.ENCODE_TFBS_by_TSS_Number[TSS_Number]

            for TFBS in Current_ENCODE_TFBS:
                # continue
                ### TFBS will be divided into 20 fragments fo get the best score value
                Fragment_size = ((TFBS.End - TFBS.Start)/20)
                ScoringPositions = [Fragment_size*i + TFBS.Start for i in range(21)]
                AdjustedScore = max([AdjustScoreAccordingToTSS_dist_ENCODE(TFBS.Score,(-1 + 2*(G.Strand == '+'))*(x - TSS)) for x in ScoringPositions])
                Replace_old_Gene_TF_feature = False
                if TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position.keys():
                    if TFBS.Start in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name].keys():
                        if AdjustedScore <= Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start]:
                            continue
                        else: Replace_old_Gene_TF_feature = True

                if not GeneName in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name.keys():
                    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName] = {}
                # print('------------------')
                # print(TFBS.Start)
                # print(TFBS.End)
                # print(TFBS.TF_Name)
                # print(TFBS.Score)
                # print(AdjustedScore)
                # print(TSS)
                # print('------------------')
                if (not TFBS.TF_Name in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName]) or Replace_old_Gene_TF_feature:  ### if this is first TFBS for this TF in the current gene
                    if Replace_old_Gene_TF_feature and not TFBS.TF_Name in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName]:
                        print('[10] Encode.py: Something strange.... exit and debug!!!')
                        exit()
                    Feat = TENCODE_ChIPSeq_Gene_TF_Features()
                    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName][TFBS.TF_Name] = Feat
                    Feat.AdjustedScore = AdjustedScore
                    Feat.Score = TFBS.Score
                    Feat.ExperimentsCount = TFBS.ExperimentsCount
                    if G.SegDescriptions_by_TSS_Number[TSS_Number] != []:
                        SegDescriptions_for_current_area = GetSegDescriptionsForSpecificArea(G.SegDescriptions_by_TSS_Number[TSS_Number],G.Chr,TFBS.Start,TFBS.End,refine = True)
                        SegDes_frequency_by_SegDes_type = {}
                        for SegDes in SegDescriptions_for_current_area:
                            for Segment in SegDes.SegListsByChr[G.Chr] if G.Chr in SegDes.SegListsByChr.keys() else []:
                                if Segment.Type in SegDes_frequency_by_SegDes_type.keys(): SegDes_frequency_by_SegDes_type[Segment.Type] += 1
                                else: SegDes_frequency_by_SegDes_type[Segment.Type] = 1

                        Feat.Genome_States_Count_by_Type = SegDes_frequency_by_SegDes_type
                    if not TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position: Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name] = {}
                    Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start] = AdjustedScore

                else:      ## if this is second, third etc. TFBS for this TF in the current gene
                    Feat = ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName][TFBS.TF_Name]
                    Feat.AdjustedScore = (Feat.AdjustedScore**2 + AdjustedScore**2)**0.5
                    Feat.Score = TFBS.Score = (Feat.Score**2 + TFBS.Score**2)**0.5
                    Feat.ExperimentsCount += TFBS.ExperimentsCount
                    if G.SegDescriptions_by_TSS_Number[TSS_Number] != []:
                        SegDescriptions_for_current_area = GetSegDescriptionsForSpecificArea(G.SegDescriptions_by_TSS_Number[TSS_Number],G.Chr,TFBS.Start,TFBS.End,refine = True)
                        for SegDes in SegDescriptions_for_current_area:
                            for Segment in SegDes.SegListsByChr[G.Chr] if G.Chr in SegDes.SegListsByChr.keys() else []:
                                if Segment.Type in Feat.Genome_States_Count_by_Type.keys(): Feat.Genome_States_Count_by_Type[Segment.Type] += 1
                                else: Feat.Genome_States_Count_by_Type[Segment.Type] = 1
                    if not TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position: Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name] = {}
                    Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start] = AdjustedScore

        GeneN += 1

    # Feat = ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name['COL16A1']['BCL11A']
    # print(Feat.Score)
    # print(Feat.AdjustedScore)
    # print(Feat.ExperimentsCount)
    # print(Feat.Genome_States_Count_by_Type)
    # exit()

    # TestList = Get_ENCODE_ChIPSeq_TFBS_for_specified_interval(ENCODE_ChIPSeq_TFBS_by_Chr,'1',32173317,32173397,Flashback_distance = 50000)
    # for x in TestList:
    #     print(x.TF_Name)
    #     print(x.Chr)
    #     print(x.Start)
    #     print(x.End)
    #     print(x.Score)
    #     print(x.ExperimentsCount)
    #     print('\n\n\n')

    return ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name


def Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name(GenesByName,ENCODE_ChIPSeq_TFBS_by_Chr,LookupDistance_upstream=25000,LookupDistance_downstream=4000,Flashback_distance=15000,SegDescriptions=[]):
    print('Performing ENCODE ChIP-Seq TFBS assignment to genes and calculating scores')
    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name = {}
    GeneN = 0
    for GeneName in sorted(GenesByName.keys()):
         # AL672183.2  problem gene
        G = GenesByName[GeneName]
        sys.stdout.write('\rProc.gene %s, %d of %d...'%(G.Name,GeneN+1,len(GenesByName)))

        TSSs = []
        if G.Strand == '+':
            TSSs = sorted(set([T.ExonStarts[0] for T in G.Transcripts]))
        elif G.Strand == '-':
            TSSs = sorted(set([T.ExonEnds[-1] for T in G.Transcripts]))
        else:
            print('[8] Smth strange. Gene %s. Maybe incorrect GTF file.'%(G.Name))
        # MinCoord = None
        # MaxCoord = None
        Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position = {}
        for TSS in TSSs:
            if G.Strand == '+':
                Current_LookupDistance_upstream = min(LookupDistance_upstream,TSS)
                Current_LookupDistance_downstream =  LookupDistance_downstream
            else:
                Current_LookupDistance_upstream = min(LookupDistance_downstream,TSS)
                Current_LookupDistance_downstream = LookupDistance_upstream
            MinCoord = TSS - Current_LookupDistance_upstream
            MaxCoord = TSS + Current_LookupDistance_downstream
            Current_ENCODE_TFBS = Get_ENCODE_ChIPSeq_TFBS_for_specified_interval(ENCODE_ChIPSeq_TFBS_by_Chr,G.Chr,MinCoord,MaxCoord,Flashback_distance)


            for TFBS in Current_ENCODE_TFBS:
                # continue
                ### TFBS will be divided into 20 fragments fo get the best score value
                Fragment_size = ((TFBS.End - TFBS.Start)/20)
                ScoringPositions = [Fragment_size*i + TFBS.Start for i in range(21)]
                AdjustedScore = max([AdjustScoreAccordingToTSS_dist_ENCODE(TFBS.Score,(-1 + 2*(G.Strand == '+'))*(x - TSS)) for x in ScoringPositions])
                Replace_old_Gene_TF_feature = False
                if TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position.keys():
                    if TFBS.Start in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name].keys():
                        if AdjustedScore <= Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start]:
                            continue
                        else: Replace_old_Gene_TF_feature = True

                if not GeneName in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name.keys():
                    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName] = {}
                # print('------------------')
                # print(TFBS.Start)
                # print(TFBS.End)
                # print(TFBS.TF_Name)
                # print(TFBS.Score)
                # print(AdjustedScore)
                # print(TSS)
                # print('------------------')
                if (not TFBS.TF_Name in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName]) or Replace_old_Gene_TF_feature:  ### if this is first TFBS for this TF in the current gene
                    if Replace_old_Gene_TF_feature and not TFBS.TF_Name in ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName]:
                        print('[10] Encode.py: Something strange.... exit and debug!!!')
                        exit()
                    Feat = TENCODE_ChIPSeq_Gene_TF_Features()
                    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName][TFBS.TF_Name] = Feat
                    Feat.AdjustedScore = AdjustedScore
                    Feat.Score = TFBS.Score
                    Feat.ExperimentsCount = TFBS.ExperimentsCount
                    if SegDescriptions != []:
                        SegDescriptions_for_current_area = GetSegDescriptionsForSpecificArea(SegDescriptions,G.Chr,TFBS.Start,TFBS.End,refine = True)
                        SegDes_frequency_by_SegDes_type = {}
                        for SegDes in SegDescriptions_for_current_area:
                            for Segment in SegDes.SegListsByChr[G.Chr] if G.Chr in SegDes.SegListsByChr.keys() else []:
                                if Segment.Type in SegDes_frequency_by_SegDes_type.keys(): SegDes_frequency_by_SegDes_type[Segment.Type] += 1
                                else: SegDes_frequency_by_SegDes_type[Segment.Type] = 1

                        Feat.Genome_States_Count_by_Type = SegDes_frequency_by_SegDes_type
                    if not TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position: Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name] = {}
                    Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start] = AdjustedScore

                else:      ## if this is second, third etc. TFBS for this TF in the current gene
                    Feat = ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName][TFBS.TF_Name]
                    Feat.AdjustedScore = (Feat.AdjustedScore**2 + AdjustedScore**2)**0.5
                    Feat.Score = TFBS.Score = (Feat.Score**2 + TFBS.Score**2)**0.5
                    Feat.ExperimentsCount += TFBS.ExperimentsCount
                    if SegDescriptions != []:
                        SegDescriptions_for_current_area = GetSegDescriptionsForSpecificArea(SegDescriptions,G.Chr,TFBS.Start,TFBS.End,refine = True)
                        for SegDes in SegDescriptions_for_current_area:
                            for Segment in SegDes.SegListsByChr[G.Chr] if G.Chr in SegDes.SegListsByChr.keys() else []:
                                if Segment.Type in Feat.Genome_States_Count_by_Type.keys(): Feat.Genome_States_Count_by_Type[Segment.Type] += 1
                                else: Feat.Genome_States_Count_by_Type[Segment.Type] = 1
                    if not TFBS.TF_Name in Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position: Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name] = {}
                    Adjusted_scores_of_submitted_TFBS_by_TF_Name_by_position[TFBS.TF_Name][TFBS.Start] = AdjustedScore

        GeneN += 1

    # Feat = ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name['COL16A1']['BCL11A']
    # print(Feat.Score)
    # print(Feat.AdjustedScore)
    # print(Feat.ExperimentsCount)
    # print(Feat.Genome_States_Count_by_Type)
    # exit()

    # TestList = Get_ENCODE_ChIPSeq_TFBS_for_specified_interval(ENCODE_ChIPSeq_TFBS_by_Chr,'1',32173317,32173397,Flashback_distance = 50000)
    # for x in TestList:
    #     print(x.TF_Name)
    #     print(x.Chr)
    #     print(x.Start)
    #     print(x.End)
    #     print(x.Score)
    #     print(x.ExperimentsCount)
    #     print('\n\n\n')

    return ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name


def Load_ENCODE_ChIPSeq_V3_data(Encode_BED_FileName,GTF_FileName,LookupDistance_upstream=25000,LookupDistance_downstream=4000,Flashback_distance=15000,SegDescriptions=[]):
    GenesByName,GenesByID,TranscriptsByID,ProteinsByID = Load_GTF_for_ENCODE(GTF_FileName)
    ENCODE_ChIPSeq_TFBS_by_Chr = {}
    for cells in ReadTableFromTXT_Autoformat(Encode_BED_FileName):
        Chr = cells[1][3:]
        Start = cells[2]
        End = cells[3]
        TF_Name = cells[4]
        Score = cells[5]
        ExperimentsCount = cells[6]
        if type(Start) != int or type(End) != int or type(Score) != int or type(ExperimentsCount) != int:
            print('Incorrect ENCODE ChIPSeq file, line %s'%cells)
            exit()
        if not Chr in ENCODE_ChIPSeq_TFBS_by_Chr.keys(): ENCODE_ChIPSeq_TFBS_by_Chr[Chr] = []
        ENCODE_ChIPSeq_TFBS_by_Chr[Chr].append(TENCODE_ChIPSeq_TFBS(TF_Name,Chr,Start,End,Score,ExperimentsCount))

    for Chr in sorted(ENCODE_ChIPSeq_TFBS_by_Chr.keys()):
        ENCODE_ChIPSeq_TFBS_by_Chr[Chr].sort(key = (lambda x: x.Start))

    SupplyGenesWith_segment_and_TFBS_info_from_ENCODE(GenesByName,ENCODE_ChIPSeq_TFBS_by_Chr,Flashback_distance,LookupDistance_upstream,LookupDistance_downstream,SegDescriptions)

    # ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name = Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name_Multiprocess(GenesByName)
    # return ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name


    ThreadsCount = min(32,min(cpu_count(),max(4,cpu_count()-4)))
    ThreadsCount = CorrectThreadCountAccordingToAvailableMem(ThreadsCount)
    pool = Pool(processes = ThreadsCount)

    StartupArray = []
    for T in range(ThreadsCount):
        Pars = Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name_Multiprocess_StartupPars(T,{})
        StartupArray.append(Pars)
    GeneN = 0
    for GeneName in sorted(GenesByName.keys()):
        ThreadN = GeneN % ThreadsCount
        StartupArray[ThreadN].GenesByName[GeneName] = GenesByName[GeneName]
        GeneN += 1

    print('Running ENCODE ChIP-Seq assignment and score calculations @%d cores...'%(ThreadsCount))
    ResultsArray = pool.map(Get_ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name_Multiprocess, StartupArray)
    pool.close()
    pool.join()
    print('\rCompleted.                                    ')
    print('Integrating ENCODE ChIP-Seq assignment and score calculations from @%d threads...'%(ThreadsCount))

    ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name = {}
    for R in ResultsArray:
        for GeneName in R.keys():
            ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name[GeneName] = R[GeneName]

    print('Completed.')

    return ENCODE_ChIPSeq_Features_by_GeneName_by_TF_name


# ActivePromoter_Color = (255,255,0) #'#ff0000'
# PromoterFlanking_Color = (249,128,128) #'#f98080'
# InactivePromoter_Color = (232,157,255) # '#e89dff'
# StrongEnhancer_Color = (250,202,0) '#faca00'
ActivePromoter_Color = '#ff0000'
PromoterFlanking_Color = '#f98080'
InactivePromoter_Color = '#e89dff'
StrongEnhancer_Color = '#faca00'
WeakEnhancer_Color = '#fff29b'
Insulator_Color = '#82aafa'
TranscriptionAssociated_Color = '#01b954'
LowActivity_Color = '#a8fca8'
PolycombRepressed_Color = '#9d9d9d'
Heterochromatin_Color = '#f5f5f5'

ColorsBy_ENCODE_SegState = dict()
ColorsBy_ENCODE_SegState['TSS'] = ActivePromoter_Color
ColorsBy_ENCODE_SegState['Tss'] = ActivePromoter_Color
ColorsBy_ENCODE_SegState['DnaseD'] = ActivePromoter_Color
ColorsBy_ENCODE_SegState['PF'] = PromoterFlanking_Color
ColorsBy_ENCODE_SegState['TssF'] = PromoterFlanking_Color
ColorsBy_ENCODE_SegState['PromF'] = PromoterFlanking_Color
ColorsBy_ENCODE_SegState['PromP'] = InactivePromoter_Color
ColorsBy_ENCODE_SegState['Enh'] = StrongEnhancer_Color
ColorsBy_ENCODE_SegState['EnhF'] = StrongEnhancer_Color
ColorsBy_ENCODE_SegState['EnhPr'] = StrongEnhancer_Color
ColorsBy_ENCODE_SegState['EnhP'] = StrongEnhancer_Color
ColorsBy_ENCODE_SegState['E'] = StrongEnhancer_Color
ColorsBy_ENCODE_SegState['WE'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['EnhWF'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['EnhW'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['DNaseU'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['DNaseD'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['FaireW'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['EnhWf'] = WeakEnhancer_Color
ColorsBy_ENCODE_SegState['CTCF'] = Insulator_Color
ColorsBy_ENCODE_SegState['CtrcfO'] = Insulator_Color
ColorsBy_ENCODE_SegState['Ctcf'] = Insulator_Color
ColorsBy_ENCODE_SegState['CtcfO'] = Insulator_Color
ColorsBy_ENCODE_SegState['T'] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['Gen5\''] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['Gen3\''] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['Elon'] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['ElonW'] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['Pol2'] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['H4K20'] = TranscriptionAssociated_Color
ColorsBy_ENCODE_SegState['Low'] = LowActivity_Color
ColorsBy_ENCODE_SegState['R'] = PolycombRepressed_Color
ColorsBy_ENCODE_SegState['ReprD'] = PolycombRepressed_Color
ColorsBy_ENCODE_SegState['Repr'] = PolycombRepressed_Color
ColorsBy_ENCODE_SegState['ReprW'] = PolycombRepressed_Color
ColorsBy_ENCODE_SegState['Quies'] = Heterochromatin_Color
ColorsBy_ENCODE_SegState['Art'] = Heterochromatin_Color
ColorsBy_ENCODE_SegState['Quiesc'] = Heterochromatin_Color
