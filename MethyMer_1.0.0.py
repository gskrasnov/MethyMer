# -*- coding: utf-8 -*-

"""
Created on Wed Feb 12 19:46:49 2014

@author: Administrator
"""
from __future__ import print_function
import math
import optparse,sys
import os
import gc
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna,generic_protein,generic_nucleotide
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# from Bio import SeqUtils
from math import fabs
import pickle
import hashlib
from time import gmtime, strftime
import textwrap
import numpy as np
# from scipy.interpolate import interp1d
import scipy.ndimage
import random
import multiprocessing
# from multiprocessing import Process
from multiprocessing import Pool
from Thermal import *
from Energy import *
from Encode import *
from Exported_Classes import *
import openpyxl

# print(str(Seq('GGAAGAGGTTTTTTTGGGGTTTTCG').complement()))
# exit()

# if 'win' in sys.platform:

import PIL.Image as Image
import PIL.ImageFont as ImageFont
import PIL.ImageDraw as ImageDraw
import PIL.ImageFont as ImageFont
    # from tkinter import Tk, Canvas, Frame, BOTH, NW,NE
    # top = Tk()

#### градиент   нормальный, без серых областей переходных
HypermethylationGradientScorePoints = [(138,222,83),(221,224,0),(255,15,15)]

# Это старый градиент, неправильный
# HypermethylationGradientScorePoints = [(98,90,255),(221,224,0),(255,10,10)]

#### градиент   синий - зелёный - розовый
HypomethylationGradientScorePoints = [(33,242,255),(92,103,255)]

class TTargetRegionCoordinates():
    def __init__(self,Gene,Chr,start,fin,Compl,IsFragmentOfLargeGene = False):
        self.Gene = Gene
        self.Chr = Chr
        self.start = start
        self.fin = fin
        self.Compl = Compl
        self.IsFragmentOfLargeGene = IsFragmentOfLargeGene


def Gradient(Points,Position):
    Position = min(1,max(0,Position))
    if Position == 1:
        return Points[-1]
    FragmentSize = (1/(len(Points)-1))
    FragmentNumber = int(Position / FragmentSize)

    InFragmentPosition = (Position%FragmentSize)/FragmentSize

    RColor=int(Points[FragmentNumber][0] + (Points[FragmentNumber+1][0] - Points[FragmentNumber][0])*InFragmentPosition)
    GColor=int(Points[FragmentNumber][1] + (Points[FragmentNumber+1][1] - Points[FragmentNumber][1])*InFragmentPosition)
    BColor=int(Points[FragmentNumber][2] + (Points[FragmentNumber+1][2] - Points[FragmentNumber][2])*InFragmentPosition)

    return (RColor,GColor,BColor)

class TCpG():
    def __init__(self,Chr,Pos):
        self.Chr = Chr
        self.Pos = Pos
        self.Pooled_T_SamplesCountByFile_Methyl = {}
        self.Pooled_C_SamplesCountByFile_Methyl = {}
        self.Paired_SamplesCountByFile_Methyl = {}
        self.Pooled_T_SamplesCountByFile_RnaSeq = {}
        self.Pooled_C_SamplesCountByFile_RnaSeq = {}
        self.Paired_SamplesCountByFile_RnaSeq = {}
        self.Pooled_T_SamplesCountByFile_Both = {}
        self.Pooled_C_SamplesCountByFile_Both = {}
        self.Paired_SamplesCountByFile_Both = {}

    #  скоринг метилирования и экспрессия по каждому из приведённых файлов
        self.HyperScore = {}
        self.HypoScore = {}
        self.Bv_median_T = {}
        self.Bv_median_C = {}
        self.Rs_pooled = {}
        self.Rs_paired = {}
        self.Pv_pooled = {}
        self.Pv_paired = {}
        self.RnaSeq_LogFC_Pooled = {}
        self.RnaSeq_FDR_Pooled = {}
        self.RnaSeq_LogFC_Paired = {}
        self.RnaSeq_FDR_Paired = {}
        self.RnaSeq_RelStDev = {}



def ResizeArray(src_array,desired_length):
    initial_length = len(src_array)
    if initial_length == desired_length:
        return np.array(src_array)

    res_array = np.zeros(desired_length,dtype=float)
    part_size =  initial_length / desired_length
    for part_number in range(desired_length):
        start_coord = part_number*part_size
        end_coord = start_coord + part_size
        values = []
        weights = []

        start_fragment_weight = int(min(start_coord + 1,end_coord)) - start_coord
        if start_fragment_weight > 0:
            values.append(src_array[int(start_coord)])
            weights.append(start_fragment_weight)

        end_fragment_weight = end_coord - int(max(end_coord,start_coord+1))
        if end_fragment_weight > 0:
            values.append(src_array[min(initial_length-1, int(end_coord))])
            weights.append(end_fragment_weight)


        if start_fragment_weight <= 0 and end_fragment_weight <= 0 :
            if int(start_coord) != int(end_coord):
                print('smth strange...')
            values.append(src_array[int(start_coord)])
            weights.append(1)

        for x in range(int(start_coord)+1,int(end_coord)):
            values.append(src_array[x])
            weights.append(1)

        res_array[part_number] = sum([values[x]*weights[x] for x in range(len(values))])/sum(weights)
    return res_array


def smoothListGaussian(list,degree=5):
     window=degree*2-1
     weight=np.array([1.0]*window)
     weightGauss=[]
     for i in range(window):
         i=i-degree+1
         frac=i/float(window)
         gauss=1/(np.exp((4*(frac))**2))
         weightGauss.append(gauss)
     weight=np.array(weightGauss)*weight
     smoothed=[0.0]*(len(list)-window)
     for i in range(len(smoothed)):
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)
     return smoothed


def intWithCommas(x):
    if type(x) not in [int]:
        raise TypeError("Parameter must be an integer.")
    if x < 0:
        return '-' + intWithCommas(-x)
    result = ''
    while x >= 1000:
        x, r = divmod(x, 1000)
        result = ",%03d%s" % (r, result)
    return "%d%s" % (x, result)

def color(R,G,B):
    Rhex=hex(R)[2:]
    if len(Rhex)==1:
        Rhex='0'+Rhex
    Ghex=hex(G)[2:]
    if len(Ghex)==1:
        Ghex='0'+Ghex
    Bhex=hex(B)[2:]
    if len(Bhex)==1:
        Bhex='0'+Bhex
    return '#%s%s%s'%(Rhex,Ghex,Bhex)


def InvertLetter(Let):
    if Let=='A': return 'T'
    if Let=='T': return 'A'
    if Let=='G': return 'C'
    if Let=='C': return 'G'
    if Let=='a': return 't'
    if Let=='t': return 'a'
    if Let=='g': return 'c'
    if Let=='c': return 'g'
    return Let

def MakeComplSeq(aSeq):
    Length=len(aSeq)
    if aSeq==None: return None
    Pos=0
    NewSeq=[' ']*Length
    while Pos < Length:
        NewSeq[Pos]=InvertLetter(aSeq[Length-Pos-1])
        Pos+=1
    return ''.join(NewSeq)


def CountMismatches(SeqA,SeqB):
    Len=min(len(SeqA),len(SeqB))
    if not type(SeqA) is str:
        SeqA=str(SeqA)
    if not type(SeqB) is str:
        SeqB=str(SeqB)
    TotalMM=0
    for P in range(Len):
        if SeqA[P]!=SeqB[P]:
            TotalMM+=1
    return TotalMM

# class TSearchParameters():
#     def __init__(self):
#         self.LookupArea=0.7
#         self.MinMismatches=2
#         self.MinMismatchesInLastFive=1
        # self.MaxFalsePCR_ProductSize=3500

def hashmd5(item):
    BinVer=str(item).encode('utf_16_be')
    return hashlib.md5(BinVer).hexdigest()

class TGeneralParameters():
    def __init__(self):
        self.TotalBootstraps = 6400
        self.FalseAmpliconsToShow = 20
        self.TmCalculationMethod = 'Breslauer'
        self.SpecTestRelImpact=0.3
        self.SpecTestMaxPenalty=140
        self.TDimersdGPenaltyIncrement = -0.25
        self.FDimersdGPenaltyIncrement = 0.5
        self.CPUcount = multiprocessing.cpu_count()

    def hashmd5(self):
        CommonList=[self.TotalBootstraps,
            self.FalseAmpliconsToShow,
            self.TmCalculationMethod,
            self.SpecTestRelImpact,
            self.SpecTestMaxPenalty,
            self.TDimersdGPenaltyIncrement,
            self.FDimersdGPenaltyIncrement]
        BinVer=str(CommonList).encode('utf_16_be')
        return hashlib.md5(BinVer).hexdigest()

    def hashmd5_short(self):
        CommonList=[self.TotalBootstraps,
            self.TmCalculationMethod,
            self.TDimersdGPenaltyIncrement,
            self.FDimersdGPenaltyIncrement]
        BinVer=str(CommonList).encode('utf_16_be')
        return hashlib.md5(BinVer).hexdigest()


class TParameters():
    def __init__(self,mode):
        if mode=='SetAsDesired':
            self.MinPrimerLength=18
            self.MaxPrimerLength=26
            self.MinPrimerTm=58
            self.MaxPrimerTm=60
            self.MaxPrimerdeltaT=1
            self.MinPrimerGCContent=30
            self.MaxPrimerGCContent=70
            self.PrimerMindG=-7
            self.PairMindG=-7
            self.AllowedRepeatGNucl=3
            self.AllowedRepeatATCNucl=4
            self.AllowedGpGsitesInPrimersCount=0
            self.MinPrimerGCendContent=1
            self.MaxPrimerGCendContent=3
            self.MinAmpliconLength=170
            self.MaxAmpliconLength=300
            self.MinAmpliconGCContent=20
            self.MaxAmpliconGCContent=80
            self.MinAmpliconTm=75
            self.MaxAmpliconTm=95
            self.MinFalseAmpliconSize=2500
            self.DistanceTolerance=3
            self.MaxRepeatedRegionsOverFourLength=8

        elif mode=='SetAsLimit':
            self.MinPrimerLength=15
            self.MaxPrimerLength=33
            self.MinPrimerTm=54
            self.MaxPrimerTm=63
            self.MaxPrimerdeltaT=6
            self.MinPrimerGCContent=5
            self.MaxPrimerGCContent=90
            self.PrimerMindG=-25
            self.PairMindG=-25
            self.AllowedRepeatGNucl=5
            self.AllowedRepeatATCNucl=8
            self.AllowedGpGsitesInPrimersCount=3
            self.MinPrimerGCendContent=1
            self.MaxPrimerGCendContent=4
            self.MinAmpliconLength=120
            self.MaxAmpliconLength=340
            self.MinAmpliconGCContent=12
            self.MaxAmpliconGCContent=95
            self.MinAmpliconTm=60
            self.MaxAmpliconTm=105
            self.MinFalseAmpliconSize=50
            self.DistanceTolerance=3
            self.MaxRepeatedRegionsOverFourLength=13

        elif mode=='SetAsPenalty':
            self.MinPrimerLength=3
            self.MaxPrimerLength=2
            self.MinPrimerTm=2
            self.MaxPrimerTm=2
            self.MaxPrimerdeltaT=3
            self.MinPrimerGCContent=0.5
            self.MaxPrimerGCContent=0.5
            self.PrimerMindG=0.4
            self.PairMindG=0.4
            self.AllowedRepeatGNucl=3.3
            self.AllowedRepeatATCNucl=2
            self.AllowedGpGsitesInPrimersCount=60
            self.MinPrimerGCendContent=6
            self.MaxPrimerGCendContent=4
            self.MinAmpliconLength=0
            self.MaxAmpliconLength=0.5
            self.MinAmpliconGCContent=0
            self.MaxAmpliconGCContent=0
            self.MinAmpliconTm=0
            self.MaxAmpliconTm=0
            self.MinFalseAmpliconSize=0.045
            self.DistanceTolerance=3
            self.MaxRepeatedRegionsOverFourLength=8

        else:
            print('INCORRECT PARAMETERS CHOICE')
            exit()


    def ObjectNamesList(self):
        return ['MinPrimerLength',
        'MaxPrimerLength',
        'MinPrimerTm',
        'MaxPrimerTm',
        'MaxPrimerdeltaT',
        'MinPrimerGCContent',
        'MaxPrimerGCContent',
        'PrimerMindG',
        'PairMindG',
        'PairMindG',
        'AllowedRepeatATCNucl',
        'AllowedRepeatGNucl',
        'AllowedGpGsitesInPrimersCount',
        'MinPrimerGCendContent',
        'MaxPrimerGCendContent',
        'MinAmpliconLength',
        'MaxAmpliconLength',
        'MinAmpliconGCContent',
        'MaxAmpliconGCContent',
        'MinAmpliconTm',
        'MaxAmpliconTm'
        'MinFalseAmpliconSize',
        'DistanceTolerance']


    def hashmd5(self):
        CommonList=[self.MinPrimerLength,
                    self.MaxPrimerLength,
                    self.MinPrimerTm,
                    self.MaxPrimerTm,
                    self.MaxPrimerdeltaT,
                    self.MinPrimerGCContent,
                    self.MaxPrimerGCContent,
                    self.PrimerMindG,
                    self.PairMindG,
                    self.AllowedRepeatATCNucl,
                    self.AllowedRepeatGNucl,
                    self.AllowedGpGsitesInPrimersCount,
                    self.MinPrimerGCendContent,
                    self.MaxPrimerGCendContent,
                    self.MinAmpliconLength,
                    self.MaxAmpliconLength,
                    self.MinAmpliconGCContent,
                    self.MaxAmpliconGCContent,
                    self.MinAmpliconTm,
                    self.MaxAmpliconTm,
                    self.DistanceTolerance]
        BinVer=str(CommonList).encode('utf_16_be')
        return str(CommonList)
        return hashlib.md5(BinVer).hexdigest()


class TPosition():
    def __init__(self):
        self.Chr=''
        self.Start=0
        self.End=0
        self.MisMatches=0
        self.Gaps=0
        self.Strand=1

    def text(self):
        return '%s:%d'%(self.Chr,self.Start)

class TPrimer():
    def __init__(self):
        self.Sequence=''
        self.SequencePlus=''
        self.SequenceUnmeth=''
        self.SeqStartInSeqPlus=0
        self.SeqEndInSeqPlus=0
        self.Length=0
        self.Tm=0
        self.Start=0
        self.End=0
        self.InChrStart=0
        self.InChrEnd=0
        self.IsForward=True
        self.CpGsCount=0
        self.FalsePositivePositions={}
        self.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand={}
        self.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand={}
        self.FalsePositiveStarts_ForwardStrand_Sorted={}
        self.FalsePositiveStarts_ReverseStrand_Sorted={}
        self.Penalty=None
        self.dG=0
        self.G_Repeats=0
        self.ATC_Repeats=0
        self.GC_EndContent=0
        self.CG_content=0
        self.AvailableSecondPrimersIDs=[]
        self.CorrespondingPairsIDs=[]
        self.ID=None
        self.SpecPenalty=0
        # self.HazardPositions={}
        # self.Chr

    def CalculatePenalty(self,DesiredP,PenaltyP):
        self.Penalty = 0
        if self.Length > DesiredP.MaxPrimerLength:
            self.Penalty += (self.Length - DesiredP.MaxPrimerLength)*PenaltyP.MaxPrimerLength
        elif self.Length < DesiredP.MinPrimerLength:
            self.Penalty += (-self.Length + DesiredP.MinPrimerLength)*PenaltyP.MinPrimerLength

        if self.Tm > DesiredP.MaxPrimerTm:
            self.Penalty += (self.Tm - DesiredP.MaxPrimerTm)*PenaltyP.MaxPrimerTm
        elif self.Tm < DesiredP.MinPrimerTm:
            self.Penalty += (-self.Tm + DesiredP.MinPrimerTm)*PenaltyP.MinPrimerTm

        if self.CG_content > DesiredP.MaxPrimerGCContent:
            self.Penalty += (self.CG_content - DesiredP.MaxPrimerGCContent)*PenaltyP.MaxPrimerGCContent
        elif self.CG_content < DesiredP.MinPrimerGCContent:
            self.Penalty += (-self.CG_content + DesiredP.MinPrimerGCContent)*PenaltyP.MinPrimerGCContent

        if self.dG < DesiredP.PrimerMindG:
            self.Penalty += (DesiredP.PrimerMindG - self.dG)*PenaltyP.PrimerMindG

        if self.G_Repeats > DesiredP.AllowedRepeatGNucl:
            self.Penalty += (self.G_Repeats - DesiredP.AllowedRepeatGNucl)*PenaltyP.AllowedRepeatGNucl

        if self.ATC_Repeats > DesiredP.AllowedRepeatATCNucl:
            self.Penalty += (self.ATC_Repeats - DesiredP.AllowedRepeatATCNucl)*PenaltyP.AllowedRepeatATCNucl

        if self.CpGsCount > DesiredP.AllowedGpGsitesInPrimersCount:
            self.Penalty += (self.CpGsCount - DesiredP.AllowedGpGsitesInPrimersCount)*PenaltyP.AllowedGpGsitesInPrimersCount

        if self.GC_EndContent > DesiredP.MaxPrimerGCendContent:
            self.Penalty += (self.GC_EndContent - DesiredP.MaxPrimerGCendContent)*PenaltyP.MaxPrimerGCendContent
        elif self.GC_EndContent < DesiredP.MinPrimerGCendContent:
            self.Penalty += (-self.GC_EndContent + DesiredP.MinPrimerGCendContent)*PenaltyP.MinPrimerGCendContent

        if self.RepeatedRegionsOverFourLength > DesiredP.MaxRepeatedRegionsOverFourLength:
            self.Penalty += (self.RepeatedRegionsOverFourLength - DesiredP.MaxRepeatedRegionsOverFourLength)*PenaltyP.MaxRepeatedRegionsOverFourLength

        return self.Penalty


class TPrimerPair():
    def __init__(self):
        self.ForwardPrimer=None
        self.ReversePrimer=None
        self.dT=0
        self.dG=0
        self.AmpliconSize=0
        self.AmpliconCG_percents=-1
        self.AmpliconTm=0
        self.AmpliconSequence=''
        self.DimersInfo=None
        self.CpGsCountInPrimers=0
        self.CpGsCount=-1
        self.FalseAmplicons=[]
        self.CompleteFalseAmpliconsCount = 0
        self.BasicPenalty=None
        self.Penalty=None
        self.ID=None
        self.ObjectiveSpecTest=False
        self.NonObjectiveSpecTestInfo=''
        self.ShortID=None

    def CalculateBasicPenalty(self,DesiredP,PenaltyP,GeneralP,Slow,Compl):
        self.BasicPenalty = 0
        if self.ForwardPrimer.Penalty==None:
            self.ForwardPrimer.CalculatePenalty(DesiredP,PenaltyP)
        if self.ReversePrimer.Penalty==None:
            self.ReversePrimer.CalculatePenalty(DesiredP,PenaltyP)
        self.BasicPenalty += self.ReversePrimer.Penalty + self.ForwardPrimer.Penalty

        if Slow:
            if Compl:
                FReplacer='CA'
                RReplacer='TG'
            else:
                RReplacer='CA'
                FReplacer='TG'

            MindG=0
            FwdSet=list(set([self.ForwardPrimer.Sequence,self.ForwardPrimer.Sequence.replace('CG',FReplacer)]))
            RevSet=list(set([self.ReversePrimer.Sequence,self.ReversePrimer.Sequence.replace('CG',RReplacer)]))

            for F in FwdSet:
                for R in RevSet:
                    G_info=Energy(F,R,GeneralP.TDimersdGPenaltyIncrement,GeneralP.FDimersdGPenaltyIncrement)
                    dG=G_info.Energy
                    if dG < MindG:
                        MindG = dG
            self.dG = MindG

            if self.dG < DesiredP.PairMindG:
                self.BasicPenalty += (DesiredP.PairMindG - self.dG)*PenaltyP.PairMindG


        if self.AmpliconSize > DesiredP.MaxAmpliconLength:
            self.BasicPenalty += (self.AmpliconSize - DesiredP.MaxAmpliconLength)*PenaltyP.MaxAmpliconLength
        elif self.AmpliconSize < DesiredP.MinAmpliconLength:
            self.BasicPenalty += (-self.AmpliconSize + DesiredP.MinAmpliconLength)*PenaltyP.MinAmpliconLength

        if self.AmpliconCG_percents > DesiredP.MaxAmpliconGCContent:
            self.BasicPenalty += (self.AmpliconCG_percents - DesiredP.MaxAmpliconGCContent)*PenaltyP.MaxAmpliconGCContent
        elif self.AmpliconCG_percents < DesiredP.MinAmpliconGCContent:
            self.BasicPenalty += (-self.AmpliconCG_percents + DesiredP.MinAmpliconGCContent)*PenaltyP.MinAmpliconGCContent

        if self.AmpliconTm > DesiredP.MaxAmpliconTm:
            self.BasicPenalty += (self.AmpliconTm - DesiredP.MaxAmpliconTm)*PenaltyP.MaxAmpliconTm
        elif self.AmpliconTm < DesiredP.MinAmpliconTm:
            self.BasicPenalty += (-self.AmpliconTm + DesiredP.MinAmpliconTm)*PenaltyP.MinAmpliconTm


        return self.BasicPenalty

    def CalculateFullPenalty(self,DesiredP,PenaltyP,GeneralP,Slow,Compl):

        if self.BasicPenalty==None:
            self.CalculateBasicPenalty(DesiredP,PenaltyP,GeneralP,Slow,Compl)

        PenaltyElements=[]
        for FA in self.FalseAmplicons:
            MM_Count=FA.ForwardPrimerMismatches+FA.ReversePrimerMismatches
            if MM_Count==0:
                AddPenaltyMultiplier=130/130
            elif MM_Count==1:
                AddPenaltyMultiplier=100/130
            elif MM_Count==2:
                AddPenaltyMultiplier=85/130
            elif MM_Count==3:
                AddPenaltyMultiplier=65/130
            elif MM_Count==4:
                AddPenaltyMultiplier=47/130
            else:
                AddPenaltyMultiplier=200/(1 + pow(MM_Count,1.5)/3)/130

            if FA.Size < DesiredP.MinFalseAmpliconSize:
                PenaltyElements.append(AddPenaltyMultiplier*(DesiredP.MinFalseAmpliconSize - FA.Size))

        PenaltyPlus=(sum([x**2 for x in PenaltyElements ]))**0.5
        self.Penalty = self.BasicPenalty + min(GeneralP.SpecTestMaxPenalty,(PenaltyPlus*GeneralP.SpecTestRelImpact))
        return self.Penalty

    def CopySkeletonFrom(self,source):
        self.Penalty = source.Penalty
        self.ForwardPrimer=TPrimer()
        self.ReversePrimer=TPrimer()
        self.ForwardPrimer.Start = source.ForwardPrimer.Start
        self.ForwardPrimer.End = source.ForwardPrimer.End
        self.ForwardPrimer.CpGsCount = source.ForwardPrimer.CpGsCount
        self.ReversePrimer.Start = source.ReversePrimer.Start
        self.ReversePrimer.End = source.ReversePrimer.End
        self.ReversePrimer.CpGsCount = source.ReversePrimer.CpGsCount
        self.ShortID = source.ShortID


class TFalseAmplicon():
    def __init__(self):
        self.Chr=''
        self.Start=0
        self.End=0
        self.Size=0
        self.ForwardPrimerMismatches=0
        self.ReversePrimerMismatches=0
        self.ForwardPrimerTm=0
        self.ReversePrimerTm=0
        self.Penalty=0


class TRange():
    def __init__(self):
        self.StartIndex=None
        self.EndIndex=None

def GetRange(array,StartCoord,EndCoord):
    Range=TRange()
    Len=len(array)
    StartIndex=0
    EndIndex=Len-1
    if StartCoord > array[-1]:
        return Range
    if EndCoord < array[0]:
        return Range

    if StartCoord > array[0]:
        IStart = 0
        IEnd = Len-1
        while True:
            if IEnd - IStart <= 1:
                break
            IMiddle = int((IStart + IEnd)/2+0.5)
            if array[IMiddle] >= StartCoord:
                IEnd = IMiddle
                continue
            else:
                IStart = IMiddle
                continue

        # print(IEnd)


        while True:
            if IEnd >= 1:
                if array[IEnd-1] >= StartCoord:
                    IEnd-=1
                else:
                    break

        StartIndex = IEnd

    if EndCoord < array[-1]:
        IStart = 0
        IEnd = Len-1
        while True:
            if IEnd - IStart <= 1:
                break
            IMiddle = int((IStart + IEnd)/2+0.5)
            if array[IMiddle] >= EndCoord:
                IEnd = IMiddle
                continue
            else:
                IStart = IMiddle
                continue

        while True:
            if IStart < Len-1:
                if array[IStart+1] <= EndCoord:
                    IStart+=1
                else:
                    break

        EndIndex = IStart

    Range.StartIndex = StartIndex
    Range.EndIndex = EndIndex
    return Range

class TCompositePairsStartupParameters():
    def __init__(self):
        self.SourceSeq=None
        self.SourceCpG_coords=None
        self.Pairs=None
        self.AttemptsPerThread=None
        self.SideGap=None
        self.ThreadNumber=None
        self.Mode='Methyl'

def CompositePairs(StartupPar):
    SourceSeq=StartupPar.SourceSeq
    SourceCpG_coords=StartupPar.SourceCpG_coords
    Pairs=StartupPar.Pairs
    AttemptsPerThread=StartupPar.AttemptsPerThread
    SideGap=StartupPar.SideGap
    ThreadNumber=StartupPar.ThreadNumber

    def GetCpG_coords(coords,start,end):  # возвращает список координат GC между start и end
        ResultCoords=[]
        for C in coords:
            if C >= start and C <= end:
                ResultCoords.append(C)

        return ResultCoords

    class TUncoveredArea():
        def __init__(self):
            self.Start=0
            self.End=0


    UncoveredAreas=[]
    InitUA=TUncoveredArea()
    InitUA.Start=SideGap
    InitUA.End=len(SourceSeq)-SideGap
    UncoveredAreas.append(InitUA)

    random.seed(ThreadNumber)

    PairSequencesByAttempt=[]           # найденные последовательности пар праймеров: для каждой попытки из 20-ти.
    CumulativePenaltiesByAttempt=[]     # соответствующие штрафы для кажлой попытки из 20-ти.
    ScoresByAttempt=[]                 # соответствующие Score для каждой попытки из 20-ти
    CpGs_coveredByAttempt=[]           # чило покрытых CpG для каждой попытки из 20-ти
    CpGs_unCoveredByAttempt=[]           # чило непокрытых CpG для каждой попытки из 20-ти

    UA=UncoveredAreas[0]

    print('Thread %d started. Total %d bootstraps to proceed'%(ThreadNumber,AttemptsPerThread))
    for CurrentAttempt in range(AttemptsPerThread):
        sys.stdout.write('\rThread %d... Attempt %d of %d...'%(ThreadNumber,CurrentAttempt+1,AttemptsPerThread))
        # print('\n\n###############################################\n###############################################')
        # print('Bootstrap %d\n'%(CurrentAttempt+1))
        BootstrapPair = Pairs[random.randint(0,len(Pairs)-1)]
        CumulativePenalty=0                                     # общий штраф для данной попытки
        CpGs_coveredCount=0                                     # число покрытых CpG динуклеотидов
        CurrentCpG_coords = list([x for x in SourceCpG_coords]) # список непокрытых CpG динуклеотидов на настоящий момент


        # идём вправо - подбираем наиболее подходящую следующую пару
        CurrentPairsSequencesRight = [BootstrapPair]

        CurrentPair = BootstrapPair
        while True:
            OuterEnd = CurrentPair.ReversePrimer.End
            OuterStart = CurrentPair.ForwardPrimer.Start
            InnerEnd = CurrentPair.ReversePrimer.Start
            InnerStart = CurrentPair.ForwardPrimer.End

            CpGs_covered = GetCpG_coords(CurrentCpG_coords,InnerStart,InnerEnd)
            [CurrentCpG_coords.remove(x) for x in CpGs_covered]
            CpGs_coveredCount+=len(CpGs_covered)
            CumulativePenalty += CurrentPair.Penalty
            # print('Bootstrap %d. Pair no. %d: %d-%d has been added (mv.>>>).  %d CpG upcoverage; Penalty=%.1f'
            #       %(CurrentAttempt+1,len(CurrentPairsSequencesRight),OuterStart,OuterEnd,len(CpGs_covered),CurrentPair.Penalty))


            if CurrentPair.ForwardPrimer.CpGsCount > 0:      # если праймеры не содержат CpG, то ампликоны могут и не пересекаться
                RealStartPrev = CurrentPair.ForwardPrimer.End
            else:
                RealStartPrev = CurrentPair.ForwardPrimer.Start

            if CurrentPair.ReversePrimer.CpGsCount > 0:
                RealEndPrev = CurrentPair.ReversePrimer.Start
            else:
                RealEndPrev = CurrentPair.ReversePrimer.End

            if RealEndPrev > UA.End:
                break

            BestPair=None
            BestPairScore=None
            for P in Pairs:
                if P.ForwardPrimer.CpGsCount > 0:      # если праймеры не содержат CpG, то ампликоны могут и не пересекаться
                    RealStart = P.ForwardPrimer.End
                else:
                    RealStart = P.ForwardPrimer.Start

                if P.ReversePrimer.CpGsCount > 0:
                    RealEnd = P.ReversePrimer.Start
                else:
                    RealEnd = P.ReversePrimer.End

                if RealEndPrev >= RealEnd:
                    continue

                # BypassedFragment = SourceSeq[RealEndPrev:RealStart]    # замедляет работу
                BypassedFragmentLength = max(RealStart - RealEndPrev,0)
                BypassedCpGsCount = len(GetCpG_coords(CurrentCpG_coords,RealEndPrev,RealStart))

                # AquiredFragment = SourceSeq[max(RealEndPrev,RealStart):RealEnd]
                AquiredFragmentLength = RealEnd - max(RealEndPrev,RealStart)
                AquiredCpGcCount = len(GetCpG_coords(CurrentCpG_coords,max(RealEndPrev,RealStart),RealEnd))

                if StartupPar.Mode == 'Methyl':
                    PairScore = (AquiredCpGcCount - 3*BypassedCpGsCount-4*(P.ForwardPrimer.CpGsCount+P.ReversePrimer.CpGsCount))/((20+P.Penalty)**0.8)
                else:
                    PairScore = (AquiredFragmentLength - 3*BypassedFragmentLength)/((20+P.Penalty)**0.8)

                if BestPairScore!=None:
                    if BestPairScore < PairScore:
                        BestPairScore = PairScore
                        BestPair = P
                else:
                    BestPairScore = PairScore
                    BestPair = P

            if BestPair!=None:
                CurrentPair = BestPair
                CurrentPairsSequencesRight.append(BestPair)
            else:
                break


        # идём влево - подбираем наиболее подходящую следующую пару


        CurrentPair = BootstrapPair
        CurrentPairsSequencesLeft = [BootstrapPair]
        FirstRound=True
        while True:
            OuterEnd = CurrentPair.ReversePrimer.End
            OuterStart = CurrentPair.ForwardPrimer.Start
            InnerEnd = CurrentPair.ReversePrimer.Start
            InnerStart = CurrentPair.ForwardPrimer.End
            if not FirstRound:
                CpGs_covered = GetCpG_coords(CurrentCpG_coords,InnerStart,InnerEnd)
                [CurrentCpG_coords.remove(x) for x in CpGs_covered]
                CpGs_coveredCount+=len(CpGs_covered)
                CumulativePenalty += CurrentPair.Penalty
                # print('Bootstrap %d. Pair no. %d: %d-%d has been added (mv.<<<).  %d CpG upcoverage; Penalty=%.1f'
                #       %(CurrentAttempt,len(CurrentPairsSequencesLeft)-1,OuterStart,OuterEnd,len(CpGs_covered),CurrentPair.Penalty))
            FirstRound = False


            if CurrentPair.ForwardPrimer.CpGsCount > 0:      # если праймеры не содержат CpG, то ампликоны могут и не пересекаться
                RealStartPrev = CurrentPair.ForwardPrimer.End
            else:
                RealStartPrev = CurrentPair.ForwardPrimer.Start

            if CurrentPair.ReversePrimer.CpGsCount > 0:
                RealEndPrev = CurrentPair.ReversePrimer.Start
            else:
                RealEndPrev = CurrentPair.ReversePrimer.End

            if RealStartPrev <= UA.Start:
                # print(RealStartPrev)
                # print(UA.Start)
                # print(UA.End)
                # print(SourceSeq)
                # print(InChrStart)
                # print(InChrEnd)
                break

            BestPair=None
            BestPairScore=None
            for P in Pairs:
                if P.ForwardPrimer.CpGsCount > 0:      # если праймеры не содержат CpG, то ампликоны могут и не пересекаться
                    RealStart = P.ForwardPrimer.End
                else:
                    RealStart = P.ForwardPrimer.Start

                if P.ReversePrimer.CpGsCount > 0:
                    RealEnd = P.ReversePrimer.Start
                else:
                    RealEnd = P.ReversePrimer.End

                if RealStartPrev <= RealStart:
                    continue

                # BypassedFragment = SourceSeq[RealEnd:RealStartPrev]    # замедляет работу
                BypassedFragmentLength = max(RealStartPrev - RealEnd,0)
                BypassedCpGsCount = len(GetCpG_coords(CurrentCpG_coords,RealEnd,RealStartPrev))

                # AquiredFragment = SourceSeq[RealStart:min(RealEnd,RealStartPrev)]
                AquiredFragmentLength = min(RealEnd,RealStartPrev) - RealStart
                AquiredCpGcCount = len(GetCpG_coords(CurrentCpG_coords,RealStart,min(RealEnd,RealStartPrev)))

                if StartupPar.Mode == 'Methyl':
                    PairScore = (AquiredCpGcCount - 3*BypassedCpGsCount-4*(P.ForwardPrimer.CpGsCount+P.ReversePrimer.CpGsCount))/((20+P.Penalty)**0.8)
                else:
                    PairScore = (AquiredFragmentLength - 3*BypassedFragmentLength)/((20+P.Penalty)**0.8)

                if BestPairScore!=None:
                    if BestPairScore < PairScore:
                        BestPairScore = PairScore
                        BestPair = P
                else:
                    BestPairScore = PairScore
                    BestPair = P

            if BestPair!=None:
                CurrentPair = BestPair
                CurrentPairsSequencesLeft.append(BestPair)
            else:
                break

        # пары выбраны. Подсчёт Score и прочее

        CurrentPairsSequence=[]
        for P in list(reversed(CurrentPairsSequencesLeft)):
            CurrentPairsSequence.append(P)

        for P in CurrentPairsSequencesRight[1:]:
            CurrentPairsSequence.append(P)

        Penalty = 100*(len(CurrentPairsSequence)) + CumulativePenalty  # Штраф на число пар
        Score = 1000*(CpGs_coveredCount - 3*len(CurrentCpG_coords)) / Penalty


        PairSequencesByAttempt.append(CurrentPairsSequence)
        CumulativePenaltiesByAttempt.append(CumulativePenalty)
        CpGs_coveredByAttempt.append(CpGs_coveredCount)
        CpGs_unCoveredByAttempt.append(len(CurrentCpG_coords))
        ScoresByAttempt.append(Score)

    BestAttemptN = ScoresByAttempt.index(max(ScoresByAttempt))
    CompletePairsList = PairSequencesByAttempt[BestAttemptN]
    BestScore = ScoresByAttempt[BestAttemptN]

    print('Thread %d finished. Best Score=%.1f'%(ThreadNumber,BestScore))

    return [x.ShortID for x in CompletePairsList],BestScore


def ToBool(word):
    word = word.casefold()
    if word in ['yes','y','on']: return  True
    elif word in ['no','n','off']: return  False
    print('Incorrect input "%s"'%word)
    exit(127)


def main():
    parser = optparse.OptionParser()
    parser.add_option('--mode', dest='Mode',default='MethylPrimers', help='MethylPrimers|SimplePrimers')
    parser.add_option('--genome-fasta', dest='GenomeFileName',default=None, help='')
    parser.add_option('--genome-version', dest='GenomeVersion',default='hg19', help='')
    parser.add_option('--encode-css', dest='ENCODEChromStateSegmDatabase',default='ENCODE.css.hg19.db', help='')
    parser.add_option('--strand-info', dest='GenesStrandInfoDatabase',default='GenesStrandInfo.db', help='')
    # parser.add_option('--coords-file', dest='CoordsFileName',default=None, help='')
    parser.add_option('--custom-sidegap', dest='SideGap',default='100', help='')
    parser.add_option('--create-seq-files', dest='CreateSeqFiles',default='off', help='')
    parser.add_option('--min-amplicon-overlap', dest='MinAmpliconOverlap',default='40', help='')
    parser.add_option('--force-primers-assembly', dest='ForcePrimersFile',default=None, help='')
    parser.add_option('--bypass-spec-test', dest='BypassSpecificityTest',default='off', help='')
    parser.add_option('--custom-db-dir', dest='PrimerDatabaseFolder',default='./InMethyl_DB', help='')
    # parser.add_option('--pd-like-thermo', dest='PD_likeThermo',default='on', help='')
    parser.add_option('--parameters', dest='Parameters_File',default='PrimerParameters.xlsx', help='')
    parser.add_option('--short-side-gap', dest='ShortSideGap',default='no', help='')
    parser.add_option('--max-bowtie-hits', dest='MaxBowtieGenomicHits',default='1000000', help='')
    # parser.add_option('--attempts-count-per-thread', dest='AttemptsCountPerThread',default='10', help='')
    # parser.add_option('--cpg-score-list', dest='CpG_ScoreFile',default=None, help='')
    parser.add_option('--tcga-methyl-results', dest='CrossHubTCGAMethylationResultsFiles',default=None, help='format Colon:./TCGA/colon.txt,Rectum:./TCGA/Rectum.txt....')
    parser.add_option('--tcga-crosshub-data-dir', dest='CrossHubTCGADatabaseDirectory',default='./TCGA_data', help='')
    parser.add_option('--update-tcga-cache', dest='UpdateTCGA_Cache',default='no', help='yes|no')
    parser.add_option('--only-visualization-of-tcga-data', dest='OnlyVisualizationOfTCGAData',default='no', help='yes|no')
    parser.add_option('--add-genes-with-coords-from-tcga', dest='AddGenesWithCoordsFromTCGA',default=None, help='TIMP3,RASSF1...')
    # parser.add_option('--consider-only-tcga-cpg-with-scores-higher-than', dest='ConsiderHighScoreTCGA_CpG',default=None, help='sample: 23 (max = 100)')
    parser.add_option('--merge-hyper-and-hypo-meth', dest='MergeTCGAHyperAndHypoMethylation',default='yes', help='')
    parser.add_option('--max-allowed-fragment-length-before-zoom', dest='MaxAllowedGeneLengthBeforeZoom',default='10000', help='')
    parser.add_option('--gtf-file', dest='GtfFile',default=None, help='')
    parser.add_option('--custom-results-name', dest='CustomResultsName',default=None, help='')
    opts, args  = parser.parse_args(sys.argv[1:])

    RunTime=strftime("%Y-%m-%d__UTC.%H-%M-%S", gmtime())

    ChrSeqs={}

    ParLimit=TParameters('SetAsLimit')
    ParDesired=TParameters('SetAsDesired')
    ParPenalty=TParameters('SetAsPenalty')
    GenPar = TGeneralParameters()
    AllowedValueNames = ['MinPrimerLength','MaxPrimerLength','MinPrimerTm','MaxPrimerTm','MaxPrimerdeltaT',
        'MinPrimerGCContent','MaxPrimerGCContent','PrimerMindG','PairMindG','AllowedRepeatGNucl',
        'AllowedRepeatATCNucl','AllowedGpGsitesInPrimersCount',
        'MaxPrimerGCendContent','MinPrimerGCendContent','MinAmpliconLength','MaxAmpliconLength',
        'MinAmpliconGCContent','MaxAmpliconGCContent','MinAmpliconTm','MaxAmpliconTm','MinFalseAmpliconSize',
        'DistanceTolerance','MaxRepeatedRegionsOverFourLength', 'FalseAmpliconsToShow',
        'TmCalculationMethod', 'SpecTestRelImpact', 'SpecTestMaxPenalty', 'TDimersdGPenaltyIncrement',
        'FDimersdGPenaltyIncrement', 'CPUcount','TotalBootstraps']

    # print(opts.PrimerDatabaseFolder)
    # exit()
    if opts.Parameters_File!=None:
        try:  wb = openpyxl.load_workbook(opts.Parameters_File)
        except:
            print('Cannot load workbook %s'%opts.Parameters_File)
            exit(127)
        if 'options' not in wb.get_sheet_names():
            print('Sheet "options" is not available in workbook %s'%(opts.Parameters_File))
            exit(127)
        sheet_options = wb['options']
        # print(wb['options'][1][0].value)
        # exit()

        for n in range(1,len(list(sheet_options))+1):
            cells = [x.value for x in sheet_options[n]]
            if cells[0] is None:  continue
            if cells[0].startswith('#'):  continue
            if cells[0] == '':  continue

            if len(cells)!=4 and len(cells)!=2:
                print('Incorrect parameters. Row %d ("%s") is discarded'%(n+1,str(cells)))
                continue

            ParName=cells[0].replace(' ','').replace('.','').replace('3','T').replace('5','F').replace('-','').replace('\'','')#.replace('GC','GC_').replace('GNucl','G_Nucl').replace('Gend','G_end').replace('GpGsites','GpG_sites').replace('ATCNucl','ATC_Nucl')

            if ParName not in AllowedValueNames and ParName not in opts.__dict__.keys():
                print('Unknown parameter "%s" at row:\n "%s".\nExiting'%(ParName,str(cells)))
                exit(127)

            Desired,Limit,Penalty=0,0,0
            two_fields=False
            if len(cells)==2:  two_fields=True
            elif cells[2] == None or cells[3] == None:  two_fields=True

            if two_fields:
                try:
                    V = int(cells[1])
                except:
                    try:
                        V = float(cells[1])
                    except:
                        V = cells[1]
                        if ParName!='TmCalculationMethod' and ParName!='CPUcount' and ParName not in opts.__dict__.keys():
                            print('Incorrect parameters. Row %s'%str(cells))
                            exit()
                if ParName=='CPUcount' and V == 'auto': continue
                if ParName in opts.__dict__.keys():  setattr(opts,ParName,V)
                else:  setattr(GenPar,ParName,V)
                continue

            for Par,N in zip([ParDesired,ParLimit,ParPenalty],[1,2,3]):
                if cells[1]!='-' and cells[1]!='':
                    try:
                        V = int(cells[N])
                    except:
                        try:
                            V = float(cells[N])
                        except:
                            print('Incorrect parameter %s (value=%s)'%(ParName,cells[N]))
                            print(cells)
                            exit()
                    try:
                        setattr(Par,ParName,V)
                    except AttributeError:
                        print('Incorrect parameter name %s'%ParName)
                        exit()

    if not os.path.exists(opts.PrimerDatabaseFolder): os.mkdir(opts.PrimerDatabaseFolder)


    if opts.CustomResultsName == None:
        if not os.path.exists('Run_%s.files'%RunTime):
            os.mkdir('Run_%s.files'%RunTime)
    else:
        if not os.path.exists('%s.files'%opts.CustomResultsName):
            os.mkdir('%s.files'%opts.CustomResultsName)

    opts.OnlyVisualizationOfTCGAData = ToBool(opts.OnlyVisualizationOfTCGAData)
    if not opts.OnlyVisualizationOfTCGAData:
        if opts.CustomResultsName == None:
            PrimerInfoFile=open('Run_%s.html'%RunTime,'w')
        else:
            PrimerInfoFile=open('%s.html'%opts.CustomResultsName,'w')
        PrimerInfoFile.write("<html>\n<body>\n<font face=\"consolas\">")
    else:
        PrimerInfoFile = None

    opts.BypassSpecificityTest = ToBool(opts.BypassSpecificityTest)
    opts.MergeTCGAHyperAndHypoMethylation = ToBool(opts.MergeTCGAHyperAndHypoMethylation)
    opts.MaxAllowedGeneLengthBeforeZoom = int(opts.MaxAllowedGeneLengthBeforeZoom)
    opts.MaxBowtieGenomicHits = int(opts.MaxBowtieGenomicHits)


    if 'methyl' in opts.Mode.casefold():
        opts.Mode = 'Methyl'
    elif 'simple' in opts.Mode.casefold():
        opts.Mode = 'Simple'
    else:
        print('Mode should be either MethylPrimers (primers for bisulfite-conv genome) or SimplePrimers (primers for simple amplification of custom genomic locus)')
        exit()

    # if opts.CoordsFileName != None:
    #     CoordsFile=open(opts.CoordsFileName,'r',buffering=10000000)
    # else:
    #     CoordsFile = None

    for opt in list(opts.__dict__):
        if getattr(opts,opt) == 'None': setattr(opts,opt,None)


    if opts.AddGenesWithCoordsFromTCGA == 'ALL':
        EngageAllTCGA_Genes = True
        AddGenesWithCoordsFromTCGA = []
    elif opts.AddGenesWithCoordsFromTCGA != None:
        AddGenesWithCoordsFromTCGA = opts.AddGenesWithCoordsFromTCGA.replace(' ','').replace(';',',').split(',')
        EngageAllTCGA_Genes = False
    else:
        AddGenesWithCoordsFromTCGA = []
        EngageAllTCGA_Genes = False

    opts.ShortSideGap = ToBool(opts.ShortSideGap)
    opts.UpdateTCGA_Cache = ToBool(opts.UpdateTCGA_Cache)

    ForcePrimers=False
    if opts.ForcePrimersFile!=None:
        ForcePrimersFile=open(opts.ForcePrimersFile,'r')
        ForcePrimers=True

    if opts.SideGap != None:
        opts.SideGap = int(opts.SideGap)

    opts.CreateSeqFiles = ToBool(opts.CreateSeqFiles)

    # if opts.ConsiderHighScoreTCGA_CpG != None:
    #     opts.ConsiderHighScoreTCGA_CpG = True
    #     try:
    #         MinTCGA_ScoreToConsider = float(opts.ConsiderHighScoreTCGA_CpG)
    #     except:
    #         print('--consider-only-tcga-cpg-with-scores-higher-than   must be a float (=min score to consider)')
    #         exit()
    # else: opts.ConsiderHighScoreTCGA_CpG = False



    # ParDesired.SpecTestRelImpact=float(opts.SpecTestRelImpact)
    # ParLimit.SpecTestRelImpact=float(opts.SpecTestRelImpact)
    # ParPenalty.SpecTestRelImpact=float(opts.SpecTestRelImpact)

    if type(GenPar.CPUcount) != int:
        try:  GenPar.CPUcount = int(GenPar.CPUcount)
        except:
            GenPar.CPUcount=multiprocessing.cpu_count()
            print('CPU count is selected as %d'%GenPar.CPUcount)

    if GenPar.TmCalculationMethod=='Breslauer':
        Thermal = TThermal(1)
        ThermalFunction=Thermal.Temper
    elif GenPar.TmCalculationMethod=='Unified':
        Thermal = TThermal(0)
        ThermalFunction=Thermal.Temper
    elif GenPar.TmCalculationMethod=='SantaLucia-mod':
        Thermal = TThermal(3)
        ThermalFunction=Thermal.Temper
    elif GenPar.TmCalculationMethod=='Sugimoto':
        Thermal = TThermal(2)
        ThermalFunction=Thermal.Temper
    elif GenPar.TmCalculationMethod=='SantaLucia':
        ThermalFunction=Tm_staluc
    else:
        print('Incorrect Tm calculation method')
        exit()

    LoadThermoParameters()



    MinAmpliconOverlap=int(opts.MinAmpliconOverlap)


    MaxScorePaired = 100
    MaxScorePooled = 200



###### Loading target region coordinates

    GeneCoordinatesList = []
    if opts.Parameters_File!=None:
        try:  wb = openpyxl.load_workbook(opts.Parameters_File)
        except:
            print('Cannot load workbook %s'%opts.Parameters_File)
            exit(127)
        if 'target regions' not in wb.get_sheet_names():
            print('Warning!! Sheet "region coordinates" is not available in workbook %s'%(opts.Parameters_File))
        else:
            sheet_coordinates = wb['target regions']
            # print(wb['options'][1][0].value)
            # exit()
            for n in range(1,len(list(sheet_coordinates))+1):
                cells = [x.value for x in sheet_coordinates[n]]
                if cells[0] is None:  continue
                if cells[0].startswith('#'):  continue
                if cells[0] == '':  continue

                if len(cells)!=5:
                    print('Incorrect parameters. Row %d ("%s") is discarded'%(n+1,str(cells)))
                    continue

                Gene=cells[0].replace('\\','_').replace('/','_').replace(':','_').replace('|','_').replace('(','_').replace(')','_')
                Chr=cells[1].replace('chr','').replace('Chr','')
                region_start=cells[2].replace(',','').replace(' ','')
                region_end=cells[3].replace(',','').replace(' ','')
                Compl = ToBool(cells[4])

                try:
                    region_start=int(region_start)
                    region_end=int(region_end)
                except ValueError:
                    print('Incorrect row %s at file %s'%(str(cells),opts.Parameters_File))
                    exit(127)
                GeneCoordinatesList.append(TTargetRegionCoordinates(Gene,Chr,region_start,region_end,Compl))
            print('Total %d target regions are specified'%(len(GeneCoordinatesList)))


####################################################################################
#######  ЗАГРУЗКА И КЭШИРОВАНИЕ ДАННЫХ О ТОМ, КАКИЕ ГЕНЫ НАХОДЯТСЯ НА ПРЯМОЙ ИЛИ ОБРАТНОЙ НИТЯХ
    GenesComplStatus = {}
    GeneComplCacheFileName = opts.GenesStrandInfoDatabase
    if os.path.exists(GeneComplCacheFileName):
        with open(GeneComplCacheFileName,'rb') as input:
            GenesComplStatus = pickle.load(input)
        print('Gene strand info is loaded from cache.')

    if opts.GtfFile != None:

        print('Loading strand info from GTF file...')
        F = open(opts.GtfFile,'r',10000000)
        while True:
            L = F.readline()
            if len(L)==0:
                break
            if len(L)<3:
                continue
            if L.startswith('#'):
                continue
            cells = L.split('\t')
            if len(cells) < 9:
                print('Incorrect GTF file')
                exit()
            if cells[2]!='exon':
                continue

            Compl = (cells[6]=='-')
            Des=(cells[-1].replace(' ','').replace(';;',';').split(';'))[:-1]
            Items={}
            for D in Des:
                Items[D.split('"')[0]]=D.split('"')[1]
            try:
                GeneName = Items['gene_name']
            except KeyError:
                continue

            GenesComplStatus[GeneName] = Compl

        print('Strand info for %d genes is found'%(len(GenesComplStatus)))
        with open(GeneComplCacheFileName,'wb') as output:
            pickle.dump(GenesComplStatus, output, pickle.HIGHEST_PROTOCOL)
        print('Gene strand info is cached')

####################################################################################
#######  Loading ENCODE chromatin state segmentation info

    SegDescriptions = []
    SegListsByCellTypeByML_TypeByChr = {}
    if os.path.exists(opts.ENCODEChromStateSegmDatabase):
        if not opts.GenomeVersion in opts.ENCODEChromStateSegmDatabase:
            print('#'*80 + '\n\n' + ' '*10 + 'Attention!! Please check consistency of the genome version (%s) and Encode state segmentation database (%s).\n'%(opts.GenomeVersion,opts.ENCODEChromStateSegmDatabase) + ' '*10 + 'A conflict may exist\n\n' + '#'*80 + '\n')
        # else:
        print('Loading ENCODE genome segments database...')
        with open(opts.ENCODEChromStateSegmDatabase,'rb') as input:
            (SegDescriptions,SegListsByCellTypeByML_TypeByChr) = pickle.load(input)
        print('Completed')
    else:
        print('#'*80 + '\n\n' + ' '*10 + 'WARNING! Encode state segmentation database is not found\n' + ' '*10 + 'This feature will be disabled\n\n' + '#'*80 + '\n')
        # for FileN in range(len(CellTypes)):
        #     CellType = CellTypes[FileN]
        #     ML_Type = ML_Types[FileN]
        #     FileName = GenomeSegmentsFileNames[FileN]
        #     if not os.path.exists(FileName):
        #         print('ENCODE genome segmentation file %s does not exist.'%FileName)
        # 
        #     G = TGenomeSegmentationDescription(FileName,CellType,ML_Type)
        #     SegListsByChr = LoadGenomeSegmentsFromFile(FileName)
        #     G.SegListsByChr = SegListsByChr
        #     SegDescriptions.append(G)
        #     try: SegListsByCellTypeByML_TypeByChr[CellType][ML_Type] = SegListsByChr
        #     except KeyError:
        #         SegListsByCellTypeByML_TypeByChr[CellType] = {}
        #         SegListsByCellTypeByML_TypeByChr[CellType][ML_Type] = SegListsByChr
        # 
        # print('caching ENCODE genome segmentation databases...')
        # with open('ENCODE.GenomeSegmentation.%s.db'%Signature,'wb') as output:
        #     pickle.dump((SegDescriptions,SegListsByCellTypeByML_TypeByChr), output, pickle.HIGHEST_PROTOCOL)
        # print('Completed.')


####################################################################################
#######  ЗАГРУЗКА ДАННЫХ МЕТИЛИРОВАНИЯ TCGA

    class TCrossHub_data():
        def __init__(self,Dataset_Abbreviation):
            self.Dataset_Abbreviation = Dataset_Abbreviation
            self.Dataset_Name = ''
            self.SamplingInfo = None
            self.RNASeq_BasicInfo_ByGeneName = None
            self.CpG_BasicInfo_ByChr_ByPos = None

    CrossHub_data_by_TCGA_code = dict()
    if opts.CrossHubTCGADatabaseDirectory is not None:
        CrossHub_DB_files = [fn for fn in os.listdir(opts.CrossHubTCGADatabaseDirectory) if fn.endswith('.db')]
        for fn in CrossHub_DB_files:
            print('Loading CrossHub TCGA database %s...'%os.path.join(opts.CrossHubTCGADatabaseDirectory,fn))
            with open(os.path.join(opts.CrossHubTCGADatabaseDirectory,fn),'rb') as input:
                SamplingInfo,RNASeq_BasicInfo_ByGeneName,CpG_BasicInfo_ByChr_ByPos = pickle.load(input)
            CrossHub_data = TCrossHub_data(SamplingInfo.Dataset_Abbreviation)
            CrossHub_data_by_TCGA_code[SamplingInfo.Dataset_Abbreviation] = CrossHub_data
            CrossHub_data.SamplingInfo = SamplingInfo
            CrossHub_data.RNASeq_BasicInfo_ByGeneName = RNASeq_BasicInfo_ByGeneName
            CrossHub_data.CpG_BasicInfo_ByChr_ByPos = CpG_BasicInfo_ByChr_ByPos
        print('Completed')

    CrossHub_TCGA_codes = sorted(CrossHub_data_by_TCGA_code)






    CpGsByChrByPos = {}
    CpGsByGeneByPos = {}
    CpGsFound = 0
    TCGA_ResultsFileDesignations = []
    if opts.CrossHubTCGAMethylationResultsFiles!=None:
        FileN = 0
        for FileComplex in opts.CrossHubTCGAMethylationResultsFiles.split(','):
            FileN += 1
            FileName = FileComplex.split(':')[1]
            FileDesignation = FileComplex.split(':')[0]
            TCGA_ResultsFileDesignations.append(FileDesignation)
            F = open(FileName,'r')
            CurrentGene = None
            CurrentChr = None
            print('Parsing TCGA results file %d of %d...'%(FileN,len(opts.CrossHubTCGAMethylationResultsFiles.split(','))))
            FileSize = os.path.getsize(FileName)
            SizeChunk = 0
            CumSize = 0
            CurrentTCGA__Pooled_SamplesCount_Methyl_T = None
            CurrentTCGA__Pooled_SamplesCount_Methyl_C = None
            CurrentTCGA__Paired_SamplesCount_Methyl = None
            while True:
                L = F.readline()
                LLen = len(L)
                CumSize+=LLen
                SizeChunk += LLen
                if SizeChunk > 90000:
                    sys.stdout.write('\rprocessed %.1f Mb of %.1f.  %d CpG annotations found'%(CumSize/1024/1024,FileSize/1024/1024,CpGsFound))
                    SizeChunk =0
                if LLen == 0:
                    break
                if LLen < 3:
                    continue
                if L.startswith('#'):
                    continue
                if not L[-1].isprintable():
                    L = L[:-1]
                if not L[-1].isprintable():
                    L = L[:-1]
                cells = L.split('\t')
                if len(cells[1]) == 0 and CurrentChr == None:
                    print('Incorrect TCGA file format')
                    exit()
                if L.startswith('$'):
                    if L.startswith('$1'):
                        CurrentTCGA__Pooled_SamplesCount_RnaSeq_T = int(cells[1])
                        CurrentTCGA__Pooled_SamplesCount_RnaSeq_C = int(cells[2])
                        CurrentTCGA__Paired_SamplesCount_RnaSeq = int(cells[3])
                        continue
                    if L.startswith('$2'):
                        CurrentTCGA__Pooled_SamplesCount_Methyl_T = int(cells[1])
                        CurrentTCGA__Pooled_SamplesCount_Methyl_C = int(cells[2])
                        CurrentTCGA__Paired_SamplesCount_Methyl = int(cells[3])
                        continue
                    if L.startswith('$3'):
                        CurrentTCGA__Pooled_SamplesCount_Both_T = int(cells[1])
                        CurrentTCGA__Pooled_SamplesCount_Both_C = int(cells[2])
                        CurrentTCGA__Paired_SamplesCount_Both = int(cells[3])
                        continue
                    continue
                if len(cells[1]) > 0:
                    CurrentGene = cells[0]
                    CurrentChr = cells[8][4:]
                    CurrentLogFC_Pooled = float(cells[1])
                    CurrentFDR_Pooled = float(cells[2])
                    try:
                        CurrentLogFC_Paired = float(cells[3])
                    except: CurrentLogFC_Paired = 0
                    try:
                        CurrentFDR_Paired = float(cells[4])
                    except:
                        CurrentFDR_Paired = 0
                    try:
                        CurrentRelStDev = float(cells[5])
                    except: CurrentRelStDev = 0
                    continue

                Position = int(cells[8])
                HyperScore = (min(float(cells[15]),MaxScorePaired)/MaxScorePaired + min(float(cells[13])**0.8 ,MaxScorePooled**0.8)/MaxScorePooled**0.8)/2
                HypoScore = (min(float(cells[18]),MaxScorePaired)/MaxScorePaired + min(float(cells[16])**0.8 ,MaxScorePooled**0.8)/MaxScorePooled**0.8)/2

                try:
                    CpGsByChrByPos[CurrentChr]
                except KeyError:
                    CpGsByChrByPos[CurrentChr] = {}

                try:
                    CpGsByChrByPos[CurrentChr][Position]
                except KeyError:
                    CpGsByChrByPos[CurrentChr][Position] = TCpG(CurrentChr,Position)

                CpGsByChrByPos[CurrentChr][Position].HyperScore[FileDesignation] = HyperScore
                CpGsByChrByPos[CurrentChr][Position].HypoScore[FileDesignation] = HypoScore
                CpGsByChrByPos[CurrentChr][Position].Bv_median_T[FileDesignation] = float(cells[11])
                CpGsByChrByPos[CurrentChr][Position].Bv_median_C[FileDesignation] = float(cells[12])
                CpGsByChrByPos[CurrentChr][Position].Rs_pooled[FileDesignation] = float(cells[20]) if cells[20] != 'nan' else 0
                CpGsByChrByPos[CurrentChr][Position].Rs_paired[FileDesignation] = float(cells[22]) if cells[22] != 'nan' else 0
                CpGsByChrByPos[CurrentChr][Position].Pv_pooled[FileDesignation] = float(cells[21]) if cells[21] != 'nan' else 1
                CpGsByChrByPos[CurrentChr][Position].Pv_paired[FileDesignation] = float(cells[23]) if cells[23] != 'nan' else 1
                CpGsByChrByPos[CurrentChr][Position].Pooled_T_SamplesCountByFile_Methyl[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_Methyl_T
                CpGsByChrByPos[CurrentChr][Position].Pooled_C_SamplesCountByFile_Methyl[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_Methyl_C
                CpGsByChrByPos[CurrentChr][Position].Paired_SamplesCountByFile_Methyl[FileDesignation] = CurrentTCGA__Paired_SamplesCount_Methyl
                CpGsByChrByPos[CurrentChr][Position].Pooled_T_SamplesCountByFile_RnaSeq[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_RnaSeq_T
                CpGsByChrByPos[CurrentChr][Position].Pooled_C_SamplesCountByFile_RnaSeq[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_RnaSeq_C
                CpGsByChrByPos[CurrentChr][Position].Paired_SamplesCountByFile_RnaSeq[FileDesignation] = CurrentTCGA__Paired_SamplesCount_RnaSeq
                CpGsByChrByPos[CurrentChr][Position].Pooled_T_SamplesCountByFile_Both[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_Both_T
                CpGsByChrByPos[CurrentChr][Position].Pooled_C_SamplesCountByFile_Both[FileDesignation] = CurrentTCGA__Pooled_SamplesCount_Both_C
                CpGsByChrByPos[CurrentChr][Position].Paired_SamplesCountByFile_Both[FileDesignation] = CurrentTCGA__Paired_SamplesCount_Both
                CpGsByChrByPos[CurrentChr][Position].RnaSeq_LogFC_Pooled[FileDesignation] = CurrentLogFC_Pooled
                CpGsByChrByPos[CurrentChr][Position].RnaSeq_FDR_Pooled[FileDesignation] = CurrentFDR_Pooled
                CpGsByChrByPos[CurrentChr][Position].RnaSeq_LogFC_Paired[FileDesignation] = CurrentLogFC_Paired
                CpGsByChrByPos[CurrentChr][Position].RnaSeq_FDR_Paired[FileDesignation] = CurrentFDR_Paired
                CpGsByChrByPos[CurrentChr][Position].RnaSeq_RelStDev[FileDesignation] = CurrentRelStDev

                try:
                    CpGsByGeneByPos[CurrentGene]
                except KeyError:
                    CpGsByGeneByPos[CurrentGene] = {}

                CpGsByGeneByPos[CurrentGene][Position] = CpGsByChrByPos[CurrentChr][Position]
                CpGsFound += 1

            sys.stdout.write('\rCompleted.  %d CpG annotations found                             \n'%(CpGsFound))
            F.close()


    #==============================================================================
    # Loading genome file
    #==============================================================================



    print('Parsing genome file....\n')

    if opts.GenomeFileName == None:
        print('Genome file must be specified!')
        exit(127)

    GenomeFile=open(opts.GenomeFileName,'r',buffering=10000000)

    for record in SeqIO.parse(GenomeFile, "fasta"):
        ChrSeqs[record.id]=record
        sys.stdout.write("\rProcessed record '%s' with length %.1f mb    " %(record.id,len(record)/1e+6))

    Keys=list(ChrSeqs.keys())
    for Chr in Keys:
        Passed=False
        try:
            int(Chr)
            Passed=True
        except:
            if Chr=='X' or Chr=='MT' or Chr=='Y':
                Passed=True
        if not Passed:
            ChrSeqs.pop(Chr)


    sys.stdout.write("\rCompleted: 100%%. Processed %d entries              "%len(ChrSeqs))


########################################################################################################
########################################################################################################
#### БСУЛЬФИТНАЯ КОНВЕРСИЯ ГЕНОМА И СОЗДАНИЕ БАЗЫ BOWTIE
########################################################################################################
########################################################################################################

    BaseFileName=os.path.join('MethBowtieDB',os.path.split(opts.GenomeFileName)[1])
    BaseFileName=BaseFileName[:BaseFileName.rfind('.')]
    # print('\n'+BaseFileName)
    if not os.path.exists('MethBowtieDB') and not opts.OnlyVisualizationOfTCGAData:
        os.mkdir('MethBowtieDB')

    BowtieIndexFileNameSample=BaseFileName+'.fwd.1.ebwt'

    if not os.path.exists(BowtieIndexFileNameSample) and not opts.OnlyVisualizationOfTCGAData:
        if not os.path.exists(BaseFileName+'.fwd.fa'):

            ForwardGenomeFile=open(BaseFileName+'.fwd.fa','w')
            BisForwardGenomeFile=open(BaseFileName+'.bis.fwd.fa','w')
            BisForwardCpGGenomeFile=open(BaseFileName+'.bis.fwd.cpg.fa','w')
            BisReverseGenomeFile=open(BaseFileName+'.bis.rew.fa','w')
            BisReverseCpGGenomeFile=open(BaseFileName+'.bis.rew.cpg.fa','w')

            for Chr in list(sorted(ChrSeqs.keys())):
                Normal=ChrSeqs[Chr].seq
                BisForward=Seq(str(Normal).replace('C','T'))
                BisForwardCpG=Seq(str(Normal).replace('CG','PG').replace('C','T').replace('PG','CG'))
                BisReverse=Seq(str(Normal).replace('G','A'))
                BisReverseCpG=Seq(str(Normal).replace('CG','CP').replace('G','A').replace('CP','CG'))
                sys.stdout.write('\rWriting Chr.%s to methyl-genome....              '%Chr)

                MySeq=SeqRecord(Normal,id='>%s.Forward'%Chr)
                SeqIO.write(MySeq,ForwardGenomeFile,"fasta")
                MySeq=SeqRecord(BisForward,id='>%s.BisForward'%Chr)
                SeqIO.write(MySeq,BisForwardGenomeFile,"fasta")
                MySeq=SeqRecord(BisForwardCpG,id='>%s.BisForwardCpG'%Chr)
                SeqIO.write(MySeq,BisForwardCpGGenomeFile,"fasta")
                MySeq=SeqRecord(BisReverse,id='>%s.BisReverse'%Chr)
                SeqIO.write(MySeq,BisReverseGenomeFile,"fasta")
                MySeq=SeqRecord(BisReverseCpG,id='>%s.BisReverseCpG'%Chr)
                # print(str(BisReverseCpG))
                SeqIO.write(MySeq,BisReverseCpGGenomeFile,"fasta")
                # exit()

            ForwardGenomeFile.close()
            BisForwardGenomeFile.close()
            BisForwardCpGGenomeFile.close()
            BisReverseGenomeFile.close()
            BisReverseCpGGenomeFile.close()

            sys.stdout.write('\rCompleted                   \n')

            os.system('bowtie-build %s.fwd.fa %s.fwd'%(BaseFileName,BaseFileName))
            os.system('bowtie-build %s.bis.fwd.fa %s.bis.fwd'%(BaseFileName,BaseFileName))
            os.system('bowtie-build %s.bis.fwd.cpg.fa %s.bis.fwd.cpg'%(BaseFileName,BaseFileName))
            os.system('bowtie-build %s.bis.rew.fa %s.bis.rew'%(BaseFileName,BaseFileName))
            os.system('bowtie-build %s.bis.rew.cpg.fa %s.bis.rew.cpg'%(BaseFileName,BaseFileName))


    if ForcePrimers:
        ForcedPrimersList=[]
        for line in ForcePrimersFile.readlines():
            if len(line)<5:
                continue
            if line[-1].isprintable()==False:
                line=line[:-1]
            if line[-1].isprintable()==False:
                line=line[:-1]
            ForcedPrimersList.append(line)

        print('Total %d forced primers'%(len(ForcedPrimersList)))


####  Adding TCGA genes

    if EngageAllTCGA_Genes:
        AddGenesWithCoordsFromTCGA = sorted(CpGsByGeneByPos.keys())
    for G in AddGenesWithCoordsFromTCGA:
        try:
            CpGsForCurrentGene = list(CpGsByGeneByPos[G].values())
        except:
            print('######  WARNING ###########################')
            print('##   Gene %s is not found within TCGA data. Bypassed.'%G)
            print('###########################################')
            continue

        CpGsForCurrentGene.sort(key=(lambda x: x.Pos))
        CpGsCoordinatesLists = [[x.Pos for x in CpGsForCurrentGene]]

        # opts.MaxAllowedGeneLengthBeforeZoom = 10000

        while True:
            AllListsAreGood = True
            NewCpGsCoordinatesLists = []
            for L in CpGsCoordinatesLists:
                if L[-1] - L[0] > opts.MaxAllowedGeneLengthBeforeZoom:
                    Mid = (L[-1] + L[0])/2
                    ListOne = []
                    ListTwo = []
                    for C in L:
                        if C < Mid:
                            ListOne.append(C)
                        else:
                            ListTwo.append(C)
                    NewCpGsCoordinatesLists.append(ListOne)
                    NewCpGsCoordinatesLists.append(ListTwo)
                    AllListsAreGood = False
                else:
                    NewCpGsCoordinatesLists.append(L)
            CpGsCoordinatesLists = NewCpGsCoordinatesLists
            if AllListsAreGood:
                break

        if len(CpGsCoordinatesLists) > 1:
            CpGsCoordinatesLists = [[x.Pos for x in CpGsForCurrentGene]] + CpGsCoordinatesLists

        FragmentNumber = 0
        for N in range(len(CpGsCoordinatesLists)):
            L = CpGsCoordinatesLists[N]
            start = max(L[0]-100,0)
            fin = min(L[-1]+100,len(ChrSeqs[CpGsForCurrentGene[0].Chr]))
            Chr = CpGsForCurrentGene[0].Chr
            try:
                Compl =  GenesComplStatus[G]
            except KeyError:
                print('\n######  WARNING ###########################')
                print('##   Gene %s strand info is not found. Assumed as positive.'%G)
                print('###########################################')
                Compl = False

            if N == 0:
                GeneCoordinatesList.append(TTargetRegionCoordinates(G+' (complete TCGA CpG list)',Chr,start,fin,Compl))
            else:
                GeneCoordinatesList.append(TTargetRegionCoordinates(G + ' (fragment %d)'%N,Chr,start,fin,Compl,IsFragmentOfLargeGene = True))


    CurrentGeneID=0
    for GeneCoordinates in GeneCoordinatesList:
        Gene = GeneCoordinates.Gene
        Chr = GeneCoordinates.Chr
        start = GeneCoordinates.start
        fin = GeneCoordinates.fin
        Compl = GeneCoordinates.Compl

        if opts.SideGap != None:
            SideGap = opts.SideGap
        else:
            SideGap = min([start-1,len(ChrSeqs[Chr])-1-fin,int((ParLimit.MaxAmpliconLength+ParLimit.MinAmpliconLength)/2)])
            if opts.ShortSideGap:
                InitFragmentLength = fin-start
                SideGap = int(SideGap/2 + SideGap/2 * min(2000, max(InitFragmentLength - 500, 0)) / 2000)

        print('Gene %s (%d out of %d): offside sequence length was selected as %d'%(Gene,CurrentGeneID+1,len(GeneCoordinatesList),SideGap))
        start=start-SideGap
        fin=fin+SideGap

        CpG_positions=[]
        CurrentSeq=str(ChrSeqs[Chr][start:fin].seq)

        if Compl:
            print('\nGene seq %s is selected complement\n'%Gene)


        Pos=-1
        while True:
            Pos=CurrentSeq.find('CG',Pos+1)
            if Pos>0:
                CpG_positions.append(Pos+1)
            else:
                break

        if Compl:
            CurrentSeqBis = CurrentSeq.replace('G','A')
            CurrentSeqBisMethCpG = CurrentSeq.replace('CG','CP').replace('G','A').replace('CP','CG')
        else:
            CurrentSeqBis=CurrentSeq.replace('C','T')
            CurrentSeqBisMethCpG=CurrentSeq.replace('CG','PG').replace('C','T').replace('PG','CG')


        if opts.CreateSeqFiles and opts.Mode == 'Methyl':
            output=open('Gene_%s__src_Chr.%s___%d-%d.txt'%(Gene,Chr,start,fin),'w')
            output.write("LOCUS       Chr.%s, %d-%d\r\n"%(Chr,start,fin))
            output.write('FEATURES\r\n')
            output.write('     gene            1..%d\r\n'%(fin-start))
            output.write('     mRNA            join(1..%d)\r\n'%(fin-start))
            for P in CpG_positions:
                output.write('     variation %d    frequency 0.2\r\n'%P)
            output.write('ORIGIN\r\n1 ')
            output.write(CurrentSeq)
            output.write('\r\n//\r\n')
            output.close()

            output=open('Gene_%s__src_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin),'w')
            MySeq=SeqRecord(Seq(CurrentSeq),id='Gene_%s__src_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin))
            SeqIO.write(MySeq,output,"fasta")
            output.close()


            output=open('Gene_%s__bis_cpg_Chr.%s___%d-%d.txt'%(Gene,Chr,start,fin),'w')
            output.write("LOCUS       Chr.%s, %d-%d (bisulfite converted)\r\n"%(Chr,start,fin))
            output.write('FEATURES\r\n')
            output.write('     gene            1..%d\r\n'%(fin-start))
            output.write('     mRNA            join(1..%d)\r\n'%(fin-start))
            for P in CpG_positions:
                output.write('     variation %d    frequency 0.2\r\n'%P)
            output.write('ORIGIN\r\n1 ')
            output.write(CurrentSeqBisMethCpG)
            output.write('\r\n//\r\n')
            output.close()

            output=open('Gene_%s__bis_cpg_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin),'w')
            MySeq=SeqRecord(Seq(CurrentSeqBisMethCpG),id='Gene_%s__src_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin))
            SeqIO.write(MySeq,output,"fasta")
            output.close()



            output=open('Gene_%s__bis_all_Chr.%s___%d-%d.txt'%(Gene,Chr,start,fin),'w')
            output.write("LOCUS       Chr.%s, %d-%d (bisulfite converted)\r\n"%(Chr,start,fin))
            output.write('FEATURES\r\n')
            output.write('     gene            1..%d\r\n'%(fin-start))
            output.write('     mRNA            join(1..%d)\r\n'%(fin-start))
            for P in CpG_positions:
                output.write('     variation %d    frequency 0.2\r\n'%P)
            output.write('ORIGIN\r\n1 ')
            output.write(CurrentSeqBis)
            output.write('\r\n//\r\n')


            output.close()

            output=open('Gene_%s__bis_all_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin),'w')
            MySeq=SeqRecord(Seq(CurrentSeqBis),id='Gene_%s__src_Chr.%s___%d-%d.fas'%(Gene,Chr,start,fin))
            SeqIO.write(MySeq,output,"fasta")
            output.close()


        if opts.Mode == 'Methyl':
            SourceSeq = CurrentSeqBisMethCpG
            SourceSeqCompl=str(Seq(SourceSeq).reverse_complement())
        elif opts.Mode == 'Simple':
            SourceSeq = CurrentSeq
            SourceSeqCompl=str(Seq(SourceSeq).reverse_complement())

        SeqLen=len(SourceSeq)
        ForwardPrimers=[]
        ReversePrimers=[]
        Pairs=[]
        PairsByID={}


        #############################################################################################################################
        # Загрузка данных из хэша, если они там есть    ########################################################################################
        ####################################################################################################################################

        SeqHash=hashlib.md5(  CurrentSeqBisMethCpG.encode('utf_16_be')  ).hexdigest()
        ParametersHashListStr=str([ParLimit.hashmd5(),ParPenalty.hashmd5(),ParDesired.hashmd5(),GenPar.hashmd5()]).encode('utf_16_be')
        ParametersHash = hashlib.md5(ParametersHashListStr).hexdigest()
        # print(str([ParLimit.hashmd5(),ParPenalty.hashmd5(),ParDesired.hashmd5(),GenPar.hashmd5()]))
        # print(ParametersHash)
        # print(ParametersHashListStr)
        # exit()
        ParametersHashListStr_Short=str([ParLimit.hashmd5(),ParPenalty.hashmd5(),ParDesired.hashmd5(),GenPar.hashmd5_short()]).encode('utf_16_be')
        ParametersHashShort = hashlib.md5(ParametersHashListStr).hexdigest()

        PickleDB_FileNameBase=os.path.join(opts.PrimerDatabaseFolder,'%s.%s.%s'%(Gene,ParametersHash,SeqHash))


        InChrStart=start
        InChrEnd=fin

        BulkForwardPrimers=[]
        BulkReversePrimers=[]
        ForwardPrimers=[]
        ReversePrimers=[]
        DesiredForwardPrimers=[]
        DesiredReversePrimers=[]

        DumpExists=[os.path.exists(PickleDB_FileNameBase+'.level%d.db'%x) for x in [1,2,3]]


        PrimersByID={}
        if opts.OnlyVisualizationOfTCGAData:
            SpecTestForwardPrimers = []
            SpecTestReversePrimers = []
            ForwardPrimers = []
            ReversePrimers = []
            DesiredForwardPrimers = []
            DesiredReversePrimers = []
            BulkForwardPrimers = []
            BulkReversePrimers = []
            PrimersByID ={}
            Pairs =[]
            PairsByID = {}
        elif ForcePrimers:
            print('List of the primers is gained.')
            CurrentID=0
            for FP in ForcedPrimersList:
                Len=len(SourceSeq)
                Pos = SourceSeq.find(FP)
                if Pos >= 0:
                    Primer=TPrimer()
                    Primer.Start = Pos
                    Primer.End = Pos+len(FP)-1
                    Primer.InChrStart = Pos+InChrStart
                    Primer.InChrEnd = Pos+len(FP)-1+InChrStart
                    Primer.CpGsCount=FP.count('CG')
                    Primer.Chr=Chr
                    Primer.Sequence=FP
                    Primer.Tm=ThermalFunction(FP)
                    ForwardPrimers.append(Primer)
                    PrimersByID[CurrentID]=Primer
                    CurrentID+=1

                Pos = SourceSeqCompl.find(FP)
                if Pos >= 0:
                    Primer=TPrimer()
                    Primer.Start = Len-Pos-1
                    Primer.End = Len-Pos+len(FP)
                    Primer.InChrStart = PosPrimer.Start+InChrStart
                    Primer.InChrEnd = Primer.End+InChrStart
                    Primer.CpGsCount=FP.count('CG')
                    Primer.Chr=Chr
                    Primer.Sequence=FP
                    Primer.Tm=ThermalFunction(FP)
                    ReversePrimers.append(Primer)
                    PrimersByID[CurrentID]=Primer
                    CurrentID+=1
            print('%d primers are found as FWD; %d are found as REV'%(len(ForwardPrimers),len(ReversePrimers)))


            def GetKey(item):
                return item.Start

            ForwardPrimers.sort(key=GetKey)
            ReversePrimers.sort(key=GetKey)

            Pairs=[]
            if len(ForwardPrimers)!=len(ReversePrimers):
                print('Warning! lengths of FWD and REV are unequal')
            for N in range(min(len(ForwardPrimers),len(ReversePrimers))):
                NewPair=TPrimerPair()
                F=ForwardPrimers[N]
                R=ReversePrimers[N]
                NewPair.ForwardPrimer=F
                NewPair.ReversePrimer=R
                NewPair.AmpliconSize=R.End-F.Start+1
                NewPair.AmpliconSequence=str(ChrSeqs[Chr][F.InChrStart:R.InChrEnd+1].seq)
                NewPair.AmpliconTm=ThermalFunction(NewPair.AmpliconSequence)
                NewPair.CpGsCount=NewPair.AmpliconSequence.count('CG')


        elif DumpExists[0] and (not DumpExists[1]) and (not DumpExists[2]):
            print('\nLoading from file %s...'%(PickleDB_FileNameBase+'.level1.db'))
            with open(PickleDB_FileNameBase+'.level1.db','rb',10000000) as input:
                SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID = pickle.load(input)
            print('Completed.')
        elif (not DumpExists[1]) and (not DumpExists[2]):

            #############################################################################################################################
            # Первичный поиск прямых и обратных праймеров      ########################################################################################
            ####################################################################################################################################


            FoundPrimersCount=0
            for CurrentSourceSeq,what,PrimersCatalogue,DesiredPrimersCatalogue in zip([SourceSeq,SourceSeqCompl],['forward','reverse'],[BulkForwardPrimers,BulkReversePrimers],[DesiredForwardPrimers,DesiredReversePrimers]):
                for PrimerLength in range(ParLimit.MinPrimerLength,ParLimit.MaxPrimerLength+1):
                    for Pos in range(0,SeqLen-PrimerLength):
                        PrimerSeq = CurrentSourceSeq[Pos:Pos+PrimerLength]

                        "НИЗЛЕЖАЩЕЕ НУОБХОДИМО ДЛЯ ТОГО, ЧТОБЫ ЗАХВАТИТЬ В ПРАЙМЕР КОНЦЕВОЙ CpG ТИПА ATATG и ATATGC"

                        if (what == 'forward' and not Compl) or (what == 'reverse' and Compl):
                            PrimerSeqPlus = CurrentSourceSeq[Pos:Pos+PrimerLength+1]
                            SeqStartInSeqPlus = 0
                            SeqEndInSeqPlus = PrimerLength-1

                        elif (what == 'reverse' and not Compl) or (what == 'forward' and Compl):
                            PrimerSeqPlus = CurrentSourceSeq[max(Pos-1,0):Pos+PrimerLength]
                            if Pos > 0:
                                SeqStartInSeqPlus = 1
                                SeqEndInSeqPlus = PrimerLength
                            else:
                                SeqStartInSeqPlus = 0
                                SeqEndInSeqPlus = PrimerLength

                        PrimerLength=len(PrimerSeq)
                        CG_percent=(PrimerSeq.count('C')+PrimerSeq.count('G'))/len(PrimerSeq)*100
                        if CG_percent < ParLimit.MinPrimerGCContent or CG_percent > ParLimit.MaxPrimerGCContent:
                            continue

                        GC_endContent=PrimerSeq[-5:].count('C')+PrimerSeq[-5:].count('G')
                        if GC_endContent < ParLimit.MinPrimerGCendContent or GC_endContent > ParLimit.MaxPrimerGCendContent:
                            continue
                        PrevLetter=None
                        RepeatedLettersCount=0
                        RepeatedLettersCountMax=0
                        RepeatedGLettersCountMax=0
                        RepeatedRegionsOverFourLength=0
                        for i in range(0,PrimerLength):
                            if PrevLetter==None:
                                PrevLetter=PrimerSeq[i]
                                RepeatedLettersCount=1
                                continue
                            if PrimerSeq[i]==PrevLetter:
                                RepeatedLettersCount+=1
                                if RepeatedLettersCount>RepeatedLettersCountMax:
                                    RepeatedLettersCountMax=RepeatedLettersCount
                                if PrevLetter=='G':
                                    if RepeatedLettersCount>RepeatedGLettersCountMax:
                                        RepeatedGLettersCountMax=RepeatedLettersCount
                            else:
                                if RepeatedLettersCount > 4:
                                    RepeatedRegionsOverFourLength += RepeatedLettersCount
                                RepeatedLettersCount=1
                                PrevLetter=PrimerSeq[i]

                        if RepeatedLettersCount > 4:
                            RepeatedRegionsOverFourLength += RepeatedLettersCount

                        if RepeatedLettersCountMax > ParLimit.AllowedRepeatATCNucl:
                            continue
                        if RepeatedGLettersCountMax > ParLimit.AllowedRepeatGNucl:
                            continue

                        if RepeatedRegionsOverFourLength > ParLimit.MaxRepeatedRegionsOverFourLength:
                            continue

                        PrimerTm=ThermalFunction(PrimerSeq)
                        if PrimerTm > ParLimit.MaxPrimerTm:
                            continue
                        if PrimerTm < ParLimit.MinPrimerTm:
                            continue



                        if (what == 'forward' and not Compl) or (what == 'reverse' and Compl):
                            try:
                                BisAllPrimerSeq = PrimerSeqPlus.replace('CG','TG')[SeqStartInSeqPlus:SeqEndInSeqPlus+1]
                            except:
                                print(PrimerSeqPlus)
                                print(SeqEndInSeqPlus+1)
                                exit()
                        elif (what == 'reverse' and not Compl) or (what == 'forward' and Compl):
                            BisAllPrimerSeq = PrimerSeqPlus.replace('CG','CA')[SeqStartInSeqPlus:SeqEndInSeqPlus+1]
                        else:
                            print('amthsadasdsad')
                            exit()

                        if opts.Mode == 'Methyl':
                            DimersC = MaximumEnergy([PrimerSeq,BisAllPrimerSeq],GenPar.TDimersdGPenaltyIncrement,GenPar.FDimersdGPenaltyIncrement)
                        elif opts.Mode == 'Simple':
                            DimersC = MaximumEnergy([PrimerSeq,PrimerSeq],GenPar.TDimersdGPenaltyIncrement,GenPar.FDimersdGPenaltyIncrement)

                        PrimerdG=DimersC.Energy

                        if PrimerdG<ParLimit.PrimerMindG:
                            continue

                        Primer=TPrimer()
                        Primer.Sequence=PrimerSeq
                        Primer.SequencePlus = PrimerSeqPlus
                        Primer.SeqStartInSeqPlus = SeqStartInSeqPlus
                        Primer.SequenceUnmeth = BisAllPrimerSeq
                        Primer.SeqEndInSeqPlus = SeqEndInSeqPlus
                        Primer.Length=PrimerLength
                        Primer.Tm=PrimerTm
                        if what=='forward':
                            Primer.Start=Pos
                            Primer.End=Pos+PrimerLength-1
                            Primer.InChrStart = Pos+InChrStart
                            Primer.InChrEnd = Pos+PrimerLength-1+InChrStart
                        else:
                            Primer.End=len(CurrentSourceSeq)-Pos-1
                            Primer.Start=len(CurrentSourceSeq)-Pos-PrimerLength
                            Primer.InChrEnd=len(CurrentSourceSeq)-Pos-1+InChrStart
                            Primer.InChrStart=len(CurrentSourceSeq)-Pos-PrimerLength+InChrStart


                        Primer.dG = PrimerdG
                        Primer.G_Repeats=RepeatedGLettersCountMax
                        Primer.ATC_Repeats=RepeatedLettersCountMax
                        Primer.GC_EndContent=GC_endContent
                        Primer.CG_content=CG_percent
                        Primer.RepeatedRegionsOverFourLength=RepeatedRegionsOverFourLength





                        Primer.IsForward=True
                        Primer.CpGsCount=PrimerSeqPlus.count('CG')
                        Primer.Chr=Chr
                        for ChrName in ChrSeqs.keys():
                            Primer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName]={}
                            Primer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName]={}
                        Primer.ID = FoundPrimersCount
                        PrimersByID[FoundPrimersCount]=Primer
                        FoundPrimersCount+=1
                        PrimersCatalogue.append(Primer)

                        # проверка, удовлетворяет ли праймер полностью желанным критериям

                        if ParDesired.MinPrimerLength <= Primer.Length and ParDesired.MaxPrimerLength >= Primer.Length\
                            and Primer.G_Repeats <= ParDesired.AllowedRepeatGNucl and Primer.ATC_Repeats <= ParDesired.AllowedRepeatATCNucl\
                            and Primer.Tm <= ParDesired.MaxPrimerTm and Primer.Tm >= ParDesired.MinPrimerTm\
                            and Primer.dG >= ParDesired.PrimerMindG\
                            and Primer.CG_content >= ParDesired.MinPrimerGCContent and Primer.CG_content <= ParDesired.MaxPrimerGCContent\
                            and Primer.CpGsCount <= ParDesired.AllowedGpGsitesInPrimersCount\
                            and Primer.GC_EndContent >= ParDesired.MinPrimerGCendContent and Primer.GC_EndContent <= ParDesired.MaxPrimerGCendContent:

                            DesiredPrimersCatalogue.append(Primer)

                        sys.stdout.write('\rCompleted %.1f%%. Found %d %s primers     '%(100*Pos/(SeqLen-PrimerLength)/(ParLimit.MaxPrimerLength-ParLimit.MinPrimerLength+1)+100*(PrimerLength-ParLimit.MinPrimerLength)/((ParLimit.MaxPrimerLength-ParLimit.MinPrimerLength+1)),len(PrimersCatalogue),what))

                sys.stdout.write('\rFound %d %s primers (%d are "ideal")                                    \n'%(len(PrimersCatalogue),what,len(DesiredPrimersCatalogue)))


            for Primer in ReversePrimers:
                NewEnd=SeqLen-Primer.Start-1
                NewStart=SeqLen-Primer.End-1
                Primer.Start=NewStart
                Primer.End=NewEnd

                Primer.InChrStart = Primer.Start+InChrStart
                Primer.InChrEnd = Primer.Start+PrimerLength-1+InChrStart


            def GetStart(item):
                return item.Start

            def GetPenalty(item):
                return item.Penalty

            BulkForwardPrimers.sort(key=GetStart)
            BulkReversePrimers.sort(key=GetStart)



            #############################################################################################################################
            # Выбор наилучших в каждом районе  ########################################################################################
            ####################################################################################################################################

            print('Extracting best primers by coverage (soft)...')

            MaxPrimerCoverage=4

            for FP in BulkForwardPrimers:
                FP.CalculatePenalty(ParDesired,ParPenalty)
            PrimersIDAccumulator=[]
            Pos=0
            while Pos < len(SourceSeq):
                PrimerPocket=[]
                for FP in BulkForwardPrimers:
                    if FP.Start <= Pos and FP.End >= Pos:
                        PrimerPocket.append(FP)
                if len(PrimerPocket)>0:
                    PrimerPocket.sort(key=GetPenalty)
                    PrimersIDAccumulator += list([x.ID for x in PrimerPocket[:MaxPrimerCoverage]])

                Pos+=1

            ForwardPrimers=[]
            if len(PrimersIDAccumulator)>0:
                PrimersIDs=list(set(PrimersIDAccumulator))
                for ID in PrimersIDs:
                    ForwardPrimers.append(PrimersByID[ID])
                ForwardPrimers.sort(key=GetStart)


            for FP in BulkReversePrimers:
                FP.CalculatePenalty(ParDesired,ParPenalty)
            PrimersIDAccumulator=[]
            Pos=0
            while Pos < len(SourceSeq):
                PrimerPocket=[]
                for FP in BulkReversePrimers:
                    if FP.Start <= Pos and FP.End >= Pos:
                        PrimerPocket.append(FP)
                if len(PrimerPocket)>0:
                    PrimerPocket.sort(key=GetPenalty)
                    PrimersIDAccumulator += list([x.ID for x in PrimerPocket[:MaxPrimerCoverage]])

                Pos+=1

            ReversePrimers=[]
            if len(PrimersIDAccumulator)>0:
                PrimersIDs=list(set(PrimersIDAccumulator))
                for ID in PrimersIDs:
                    ReversePrimers.append(PrimersByID[ID])
                ReversePrimers.sort(key=GetStart)


            print('Completed. %d best forward primers are chosen. %d best reverse primers are chosen'%(len(ForwardPrimers),len(ReversePrimers)))


            #############################################################################################################################
            # Поиск пар       ########################################################################################
            ####################################################################################################################################

            Pairs=[]

            FPrimerNumber=0
            PairsFound=0
            for FPrimer in ForwardPrimers:
                FPrimerNumber+=1
                sys.stdout.write('\rMatching pair for %d of %d..... Found %d pairs....'%(FPrimerNumber,len(ForwardPrimers),PairsFound))
                # if PairsFound>3000:
                #     break
                for RPrimer in ReversePrimers:
                    AmpliconSize=RPrimer.End - FPrimer.Start+1
                    if AmpliconSize > ParLimit.MaxAmpliconLength or AmpliconSize < ParLimit.MinAmpliconLength:
                        continue
                    if fabs(FPrimer.Tm - RPrimer.Tm) > ParLimit.MaxPrimerdeltaT:
                        continue

                    AmpliconSeq=SourceSeq[FPrimer.Start:RPrimer.End+1]
                    AmpliconCGpercents=(AmpliconSeq.count('C')+AmpliconSeq.count('G'))/len(AmpliconSeq)*100
                    if AmpliconCGpercents > ParLimit.MaxAmpliconGCContent or AmpliconCGpercents < ParLimit.MinAmpliconGCContent:
                        continue
                    AmpliconTm=85 #ThermalFunction(AmpliconSeq)
                    if AmpliconTm > ParLimit.MaxAmpliconTm or AmpliconTm < ParLimit.MinAmpliconTm:
                        continue


                    CurrentPrimerPair=TPrimerPair()
                    CurrentPrimerPair.ForwardPrimer=FPrimer
                    FPrimer.AvailableSecondPrimersIDs.append(RPrimer.ID)
                    FPrimer.CorrespondingPairsIDs.append(PairsFound)
                    CurrentPrimerPair.ReversePrimer=RPrimer
                    RPrimer.AvailableSecondPrimersIDs.append(FPrimer.ID)
                    RPrimer.CorrespondingPairsIDs.append(PairsFound)
                    CurrentPrimerPair.dT=fabs(FPrimer.Tm - RPrimer.Tm)
                    CurrentPrimerPair.AmpliconSequence=AmpliconSeq
                    CurrentPrimerPair.AmpliconSize=AmpliconSize
                    CurrentPrimerPair.AmpliconTm=AmpliconTm
                    CurrentPrimerPair.AmpliconCG_percents=AmpliconCGpercents
                    CurrentPrimerPair.CpGsCountInPrimers=FPrimer.CpGsCount+RPrimer.CpGsCount
                    CurrentPrimerPair.CpGsCount=AmpliconSeq.count('CG')
                    CurrentPrimerPair.ID=PairsFound
                    Pairs.append(CurrentPrimerPair)
                    PairsByID[PairsFound]=CurrentPrimerPair

                    PairsFound+=1


            ##################################################################################################################################
            # Филирование с удалением похожих пар       #####################################################################################
            ####################################################################################################################################
            print('\nPair subsampling...')
            DistanceTolerance=5
            PairsMassive = list(Pairs)     # все подряд пары

            SeedPairs=[]          # базовые пары
            SeedPairsPenalties=[]          # базовые пары

            for N in range(len(Pairs)):
                Pairs[N].CalculateBasicPenalty(ParDesired,ParPenalty,GenPar,False,Compl)


            for N in range(len(Pairs)):
                sys.stdout.write('\rProcessing pair %d of %d...'%(N+1,len(Pairs)))
                F_Start=Pairs[N].ForwardPrimer.Start
                F_End=Pairs[N].ForwardPrimer.End
                R_Start=Pairs[N].ReversePrimer.Start
                R_End=Pairs[N].ReversePrimer.End

                Found=False
                for SN in range(len(SeedPairs)):
                    SF_Start = SeedPairs[SN][0].ForwardPrimer.Start
                    SF_End = SeedPairs[SN][0].ForwardPrimer.End
                    SR_Start = SeedPairs[SN][0].ReversePrimer.Start
                    SR_End = SeedPairs[SN][0].ReversePrimer.End

                    if fabs(SF_Start -F_Start) <= DistanceTolerance\
                        and fabs(SF_End -F_End) <= DistanceTolerance\
                        and fabs(SR_Start -R_Start) <= DistanceTolerance\
                        and fabs(SR_End -R_End) <= DistanceTolerance:

                        SeedPairs[SN].append(Pairs[N])
                        SeedPairsPenalties[SN].append(Pairs[N].BasicPenalty)
                        Found=True
                        break

                if not Found:
                    SeedPairs.append([])
                    SeedPairsPenalties.append([])
                    SeedPairs[-1].append(Pairs[N])
                    SeedPairsPenalties[-1].append(Pairs[N].BasicPenalty)

            BestSeedPairs=[]
            for N in range(len(SeedPairs)):
                BestPairN = SeedPairsPenalties[N].index(min(SeedPairsPenalties[N]))
                BestSeedPairs.append(SeedPairs[N][BestPairN])


            print('\rCompleted. Efficiency is %.1f%%. %d pairs remaining.         '%((len(Pairs)-len(BestSeedPairs))/len(Pairs)*100,len(BestSeedPairs)))

            Pairs=BestSeedPairs

        #############################################################################################################################
        # Выделение  праймеров, вошедших в пары (для проверки специфичности)       ######################################################################################
        ####################################################################################################################################

            SpecTestForwardPrimers = []
            SpecTestReversePrimers = []
            PairedPrimersIDs = {}
            for P in Pairs:
                PairedPrimersIDs[P.ForwardPrimer.ID] = True
                PairedPrimersIDs[P.ReversePrimer.ID] = True

            for F in ForwardPrimers:
                try:
                    PairedPrimersIDs[F.ID]
                    SpecTestForwardPrimers.append(F)
                except KeyError:
                    None

            for F in ReversePrimers:
                try:
                    PairedPrimersIDs[F.ID]
                    SpecTestReversePrimers.append(F)
                except KeyError:
                    None

            Eff=100*(len(ForwardPrimers)-len(SpecTestForwardPrimers))/len(ForwardPrimers)
            print('Efficiency relatively to forward primers is %.1f%%: %d of %d remaining for the further analysis;'%(Eff,len(SpecTestForwardPrimers),len(ForwardPrimers)))
            Eff=100*(len(ReversePrimers)-len(SpecTestReversePrimers))/len(ReversePrimers)
            print('Efficiency relatively to reverse primers is %.1f%%: %d of %d remaining.'%(Eff,len(SpecTestReversePrimers),len(ReversePrimers)))

            print('Re-calculating penalties...\n')
            N=0
            for P in Pairs:
                N+=1
                sys.stdout.write('\rProcessing pair %d of %d'%(N,len(Pairs)))
                P.CalculateBasicPenalty(ParDesired,ParPenalty,GenPar,True,Compl)

            print('\rCompleted...')



            print('Dumping...')
            with open(PickleDB_FileNameBase+'.level1.db','wb') as output:
                pickle.dump((SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID), output, pickle.HIGHEST_PROTOCOL)
            print('Completed...')



        #############################################################################################################################
        # Проверка специфичности       ######################################################################################
        ####################################################################################################################################

        SpecPenaltyThreshold = 1
        if opts.OnlyVisualizationOfTCGAData:
            None
        elif DumpExists[1] and (not DumpExists[2]):
            print('\nLoading from file %s...'%(PickleDB_FileNameBase+'.level2.db'))
            with open(PickleDB_FileNameBase+'.level2.db','rb',10000000) as input:
                SpecPenaltyThreshold,SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID = pickle.load(input)
            print('Completed.')
        elif opts.BypassSpecificityTest and (not DumpExists[2]):
            print('Specificity test bypassed')
            for P in Pairs:
                P.CalculateFullPenalty(ParDesired,ParPenalty,GenPar,True,Compl)

            print('Dumping...')
            with open(PickleDB_FileNameBase+'.level2.db','wb') as output:
                pickle.dump((SpecPenaltyThreshold,SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID), output, pickle.HIGHEST_PROTOCOL)
            print('Completed...')
        elif not DumpExists[2]:
            if opts.Mode == 'Methyl':
                PrimersFileName='%s.%s.%s.primers.fa'%(Gene,ParametersHashShort,SeqHash)
                PrimersFile=open(PrimersFileName,'w')

                SamFileNames=['%s.%s.%s.fwd.sam'%(Gene,ParametersHashShort,SeqHash),
                              '%s.%s.%s.bis.fwd.sam'%(Gene,ParametersHashShort,SeqHash),
                              '%s.%s.%s.bis.fwd.cpg.sam'%(Gene,ParametersHashShort,SeqHash),
                              '%s.%s.%s.bis.rew.sam'%(Gene,ParametersHashShort,SeqHash),
                              '%s.%s.%s.bis.rew.cpg.sam'%(Gene,ParametersHashShort,SeqHash)]

                PrimersByID={}
                CurrentID=0
                for Type,PrimerList in zip(['forward','reverse'],[SpecTestForwardPrimers,SpecTestReversePrimers]):
                    for Primer in PrimerList:
                        PrimersByID[Primer.ID]=Primer
                        CurrentID+=1
                len(PrimersByID)

                BowtieFileExists=True
                for SFN in SamFileNames:
                    if not os.path.exists(SFN):
                        BowtieFileExists=False
                        break


                PrimerSpecPenaltiesByID={}
                print('Preparing primer files...')
                for Pr in SpecTestForwardPrimers:
                    CSequences = list(sorted(set([Pr.Sequence,Pr.SequenceUnmeth])))
                    for N in range(len(CSequences)):
                        MyRec=SeqRecord(Seq(CSequences[N]),id='%s.%d.F'%(str(Pr.ID),N),description='')
                        SeqIO.write(MyRec,PrimersFile,'fasta')
                        PrimerSpecPenaltiesByID[Pr.ID]=0

                for Pr in SpecTestReversePrimers:
                    CSequences = list(sorted(set([Pr.Sequence,Pr.SequenceUnmeth])))
                    for N in range(len(CSequences)):
                        MyRec=SeqRecord(Seq(CSequences[N]),id='%s.%d.R'%(str(Pr.ID),N),description='')
                        SeqIO.write(MyRec,PrimersFile,'fasta')
                        PrimerSpecPenaltiesByID[Pr.ID]=0

                PrimersFile.close()

                if not BowtieFileExists:
                    # if not ('win' in sys.platform):
                    print('Rinning PCR specificity test...')
                    print('running bowtie for fwd. sequence...')
                    os.system('bowtie -k %d -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.fwd %s.%s.%s.primers.fa %s.%s.%s.fwd.sam'%(opts.MaxBowtieGenomicHits,min(8,GenPar.CPUcount),ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))
                    print('running bowtie for bis.fwd. sequence...')
                    os.system('bowtie -k %d -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.bis.fwd %s.%s.%s.primers.fa %s.%s.%s.bis.fwd.sam'%(opts.MaxBowtieGenomicHits,min(8,GenPar.CPUcount),ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))
                    print('running bowtie for bis.fwd.methcpg. sequence...')
                    os.system('bowtie -k %d -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.bis.fwd.cpg %s.%s.%s.primers.fa %s.%s.%s.bis.fwd.cpg.sam'%(opts.MaxBowtieGenomicHits,min(8,GenPar.CPUcount),ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))
                    print('running bowtie for bis.rev. sequence...')
                    os.system('bowtie -k %d -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.bis.rew %s.%s.%s.primers.fa %s.%s.%s.bis.rew.sam'%(opts.MaxBowtieGenomicHits,min(8,GenPar.CPUcount),ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))
                    print('running bowtie for bis.rev.methcpg sequence...')
                    os.system('bowtie -k %d -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.bis.rew.cpg %s.%s.%s.primers.fa %s.%s.%s.bis.rew.cpg.sam'%(opts.MaxBowtieGenomicHits,min(8,GenPar.CPUcount),ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))

                else:
                    print('Using previously revealed bowtie SAM files....\n')

            elif opts.Mode == 'Simple':
                PrimersFileName='%s.%s.%s.primers.fa'%(Gene,ParametersHashShort,SeqHash)
                PrimersFile=open(PrimersFileName,'w')

                SamFileNames=['%s.%s.%s.fwd.sam'%(Gene,ParametersHashShort,SeqHash)]

                PrimersByID={}
                CurrentID=0
                for Type,PrimerList in zip(['forward','reverse'],[SpecTestForwardPrimers,SpecTestReversePrimers]):
                    for Primer in PrimerList:
                        PrimersByID[Primer.ID]=Primer
                        CurrentID+=1
                len(PrimersByID)

                BowtieFileExists=True
                for SFN in SamFileNames:
                    if not os.path.exists(SFN):
                        BowtieFileExists=False
                        break


                PrimerSpecPenaltiesByID={}
                print('Preparing primer files...')
                for Pr in SpecTestForwardPrimers:
                    CSequences = list(sorted(set([Pr.Sequence,Pr.SequenceUnmeth])))
                    for N in range(len(CSequences)):
                        MyRec=SeqRecord(Seq(CSequences[N]),id='%s.%d.F'%(str(Pr.ID),N),description='')
                        SeqIO.write(MyRec,PrimersFile,'fasta')
                        PrimerSpecPenaltiesByID[Pr.ID]=0

                for Pr in SpecTestReversePrimers:
                    CSequences = list(sorted(set([Pr.Sequence,Pr.SequenceUnmeth])))
                    for N in range(len(CSequences)):
                        MyRec=SeqRecord(Seq(CSequences[N]),id='%s.%d.R'%(str(Pr.ID),N),description='')
                        SeqIO.write(MyRec,PrimersFile,'fasta')
                        PrimerSpecPenaltiesByID[Pr.ID]=0

                PrimersFile.close()

                if not BowtieFileExists:
                    # if not ('win' in sys.platform):
                    print('Rinning PCR specificity test...')
                    os.system('bowtie --all -f -v 3 --threads %d --seedlen %d --tryhard  --sam  %s.fwd %s.%s.%s.primers.fa %s.%s.%s.fwd.sam'%(GenPar.CPUcount,ParLimit.MinPrimerLength-1,BaseFileName,Gene,ParametersHashShort,SeqHash,Gene,ParametersHashShort,SeqHash))
                else:
                    print('Using previously revealed bowtie SAM files....\n')


            SamFiles=[open(x,'r') for x in SamFileNames]

            PrintedCor=False
            PrintedNCor=False


            OverMatchedPrimers={}
            PrimersSpecPenaltyByID={}
            MaxPenaltySrc = 4*opts.MaxBowtieGenomicHits
            MaxPenalty = 4*opts.MaxBowtieGenomicHits


            SamFileN=0
            for Sam in SamFiles:
                print('\nparsing sam file %d of %d...\n'%(SamFileN+1,(5 if opts.Mode == 'Methyl' else 1)))
                SizeChunk = 0
                LineN=0
                CumSize=0
                FileSize = os.path.getsize(SamFileNames[SamFileN])
                SamFileN+=1
                while True:
                    line = Sam.readline()
                    L=len(line)
                    if L==0:
                        break
                    SizeChunk+=L
                    CumSize+=L
                    LineN+=1
                    if line.startswith('@'):
                        continue
                    if L<9:
                        continue
                    if not line[-1].isprintable:
                        line=line[:-1]
                    cells=line.split('\t')
                    if len(cells[2])<=1:
                        continue
                    ID=int(cells[0].split('.')[0])
                    if SizeChunk > 90000:
                        sys.stdout.write('\rprocessed %.1f Mb of %.1f.  %d primers are have too many hits in the genome (penalty > %.1fe+6).'%(CumSize/1024/1024,FileSize/1024/1024,len(OverMatchedPrimers),MaxPenaltySrc/1000000))
                        SizeChunk =0
                    try:
                        OverMatchedPrimers[ID]
                        continue
                    except KeyError:
                        pass

                    SubSeqNumber=int(cells[0].split('.')[1])
                    PrimerType=cells[0].split('.')[2]
                    Start=int(cells[3])-1
                    SamChr=cells[2][1:].split('.')[0]
                    CurrentPosition = TPosition()
                    CurrentPosition.Start=Start
                    if cells[1]=='16':
                        CurrentPosition.Strand=2
                    elif cells[1]=='0':
                        CurrentPosition.Strand=1
                    else:
                        CurrentPosition.Strand=0
                        print('Strange alignment for primer with ID %d'%ID)
                    End=Start+len(cells[9])-1
                    CurrentPosition.End=End
                    CurrentPosition.Chr=SamChr
                    PrimerStart = PrimersByID[ID].InChrStart
                    PrimerEnd = PrimersByID[ID].InChrEnd
                    PrimerChr = PrimersByID[ID].Chr

                    if PrimerType =='F':
                        CSequences = list(sorted(set([PrimersByID[ID].Sequence,PrimersByID[ID].SequenceUnmeth])))
                    elif PrimerType =='R':
                        CSequences = list(sorted(set([PrimersByID[ID].Sequence,PrimersByID[ID].SequenceUnmeth])))
                    else:
                        print('smth strange.....')
                        exit()

                    PrimerFound = False

                    # fwd. sequence...
                    if SamFileN==1:
                        RefSequence = str(ChrSeqs[SamChr][Start:End+1].seq)

                    # bis.fwd. sequence...
                    elif SamFileN==2:
                        RefSequence = str(ChrSeqs[SamChr][Start:End+1].seq).replace('C','T')

                    # bis.fwd.methcpg. sequence...
                    elif SamFileN==3:
                        RefSequence = str(ChrSeqs[SamChr][Start:End+2].seq).replace('CG','PG').replace('C','T').replace('PG','CG')[:-1]

                    # bis.rev. sequence...
                    elif SamFileN==4:
                        RefSequence = str(ChrSeqs[SamChr][Start:End+1].seq).replace('G','A')

                    # bis.rev.methcpg sequence...
                    elif SamFileN==5:
                        if Start > 0:
                            RefSequence = str(ChrSeqs[SamChr][Start-1:End+1].seq).replace('CG','CP').replace('G','A').replace('CP','CG')[1:]
                        else:
                            RefSequence = str(ChrSeqs[SamChr][Start:End+1].seq).replace('CG','CP').replace('G','A').replace('CP','CG')



                    CurrentMismatchesCount = min(
                        CountMismatches(
                            CSequences[SubSeqNumber],
                            RefSequence
                        ),
                        CountMismatches(
                            str(Seq(CSequences[SubSeqNumber]).reverse_complement()),
                            RefSequence
                        )
                    )

                    if CurrentMismatchesCount > 6:
                        print('Smth strange.....\n   %s\n   %s\n   %s\n'%(CSequences[SubSeqNumber],str(Seq(CSequences[SubSeqNumber]).reverse_complement()),RefSequence))
                        print('---')
                        print(PrimersByID[ID].Sequence)
                        print(PrimersByID[ID].SequenceUnmeth)
                        print(PrimersByID[ID].SequencePlus)
                        print(len(PrimersByID[ID].Sequence))
                        print(PrimersByID[ID].SeqStartInSeqPlus)
                        print(PrimersByID[ID].SeqEndInSeqPlus)
                        print(line)
                        exit(130)


                    if PrimerChr == SamChr:
                        if math.fabs(PrimerStart - Start) < 4 or math.fabs(PrimerStart - End)<4 or math.fabs(PrimerEnd - Start)<4 or math.fabs(PrimerEnd - End)<3:
                            CorrectMatch=True
                        else:
                            CorrectMatch=False
                    else:
                        CorrectMatch=False

                    if not CorrectMatch:

                        AddToDict=True
                        try:
                            if CurrentPosition.Strand==1 or CurrentPosition.Strand==0:
                                PrevMisMatchCount = PrimersByID[ID].MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[SamChr][Start]
                                if PrevMisMatchCount >= CurrentMismatchesCount:
                                    AddToDict = False
                        except KeyError:
                            pass

                        try:
                            if CurrentPosition.Strand==2 or CurrentPosition.Strand==0:
                                PrevMisMatchCount = PrimersByID[ID].MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[SamChr][Start]
                                if PrevMisMatchCount >= CurrentMismatchesCount:
                                    AddToDict = False
                        except KeyError:
                            pass

                        if AddToDict:
                            PrimersByID[ID].FalsePositivePositions[CurrentPosition.text()] = CurrentPosition
                            if CurrentPosition.Strand==1 or CurrentPosition.Strand==0:
                                PrimersByID[ID].MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[SamChr][Start]=CurrentMismatchesCount

                            if CurrentPosition.Strand==2 or CurrentPosition.Strand==0:
                                PrimersByID[ID].MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[SamChr][Start]=CurrentMismatchesCount

                            PrimerSpecPenaltiesByID[ID] += max([(4 - CurrentMismatchesCount),0])

                        try:
                            PrimersSpecPenaltyByID[ID] += (5 - CurrentMismatchesCount)**2
                        except KeyError:
                            PrimersSpecPenaltyByID[ID] = (5 - CurrentMismatchesCount)**2

                        if PrimersSpecPenaltyByID[ID] > MaxPenaltySrc:
                            OverMatchedPrimers[ID] = True


            PrimerSpecPenaltiesByID[Pr.ID]=0


            [x.close() for x in SamFiles]
            # [os.remove(x) for x in SamFiles]

            for Pr in SpecTestForwardPrimers + SpecTestReversePrimers:
                for ChrName in ChrSeqs.keys():
                    Pr.FalsePositiveStarts_ForwardStrand_Sorted[ChrName] = list(Pr.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName].keys())
                    Pr.FalsePositiveStarts_ForwardStrand_Sorted[ChrName].sort()
                    Pr.FalsePositiveStarts_ReverseStrand_Sorted[ChrName] = list(Pr.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName].keys())
                    Pr.FalsePositiveStarts_ReverseStrand_Sorted[ChrName].sort()


        #######################################№№№№№№№№№######################################################################################
        # Отметка самых неспецифичных праймеров     ########################################################################################
        ################################################№№####################################################################################

            if opts.CustomResultsName == None:
                AuxInfoFile=open('Run_%s.auxinfo.txt'%RunTime,'w')
            else:
                AuxInfoFile=open('%s.auxinfo.txt'%opts.CustomResultsName,'w')

            ForwardPrimersDistributionByFalseMatches={}
            for P in SpecTestForwardPrimers:
                try:
                    ForwardPrimersDistributionByFalseMatches[len(P.FalsePositivePositions)]+=1
                except:
                    ForwardPrimersDistributionByFalseMatches[len(P.FalsePositivePositions)]=1

            AuxInfoFile.write('Distribution of FORWARD PRIMERS by false positive matches in the genome:\n')
            for K in sorted(ForwardPrimersDistributionByFalseMatches.keys()):
                AuxInfoFile.write('%d f.m.:\t%d primers\n'%(K,ForwardPrimersDistributionByFalseMatches[K]))

            AuxInfoFile.write('\n\n')


            ReversePrimersDistributionByFalseMatches={}
            for P in SpecTestReversePrimers:
                try:
                    ReversePrimersDistributionByFalseMatches[len(P.FalsePositivePositions)]+=1
                except:
                    ReversePrimersDistributionByFalseMatches[len(P.FalsePositivePositions)]=1

            AuxInfoFile.write('Distribution of REVERSE PRIMERS by false positive matches in the genome:\n\n\n')
            for K in sorted(ReversePrimersDistributionByFalseMatches.keys()):
                AuxInfoFile.write('\n%d f.m.:\t%d primers'%(K,ReversePrimersDistributionByFalseMatches[K]))


            AuxInfoFile.write('\n\n\nDistribution of FORWARD+REVERSE PRIMERS by Specificity-Penalty:\n')
            SpecPenaltyArray=[]
            for Pr in SpecTestForwardPrimers + SpecTestReversePrimers:
                SpecPenalty=0
                for CChr in Pr.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand.keys():
                    for Pos in Pr.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[CChr].keys():
                        SpecPenalty += (5-Pr.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[CChr][Pos])**2
                Pr.SpecPenalty = SpecPenalty
                SpecPenaltyArray.append(SpecPenalty)

            for SP in sorted(SpecPenaltyArray):
                AuxInfoFile.write('%d\n'%SP)


            SpecPenaltyArray.sort()
            L = len(SpecPenaltyArray)
            ThrPos=int(L*0.9)
            print('Threshold IS:  '+str(SpecPenaltyArray[ThrPos]))
            SpecPenaltyThreshold = max(min(MaxPenalty,SpecPenaltyArray[ThrPos]),100000)

            AuxInfoFile.write('\n\n\n Penalty threshold is selected as %.1f (native %.1f)'%(SpecPenaltyThreshold,SpecPenaltyArray[ThrPos]))



        #######################################№№№№№№№№№######################################################################################
        # Анализ пар праймеров     ########################################################################################
        ################################################№№####################################################################################

            print('Analyzing primer pairs...\n')
            # ThreadsCount = multiprocessing.cpu_count()
            # TotalPairs=len(Pairs)

            def GetPenalty(item):
                return item.Penalty

            PairN=0
            for PP in Pairs:
                if len(PP.FalseAmplicons) > 100000:
                    PP.ObjectiveSpecTest = False
                    if PP.NonObjectiveSpecTestInfo=='':
                        PP.NonObjectiveSpecTestInfo = 'more than 100.000 false-positive amplicons'
                    elif not PP.NonObjectiveSpecTestInfo.endswith('100.000 false-positive amplicons'):
                        PP.NonObjectiveSpecTestInfo += '; more than 100.000 false-positive amplicons'
                    continue
                PP.ObjectiveSpecTest = True
                PairN+=1
                sys.stdout.write('\ranalyzing pair %d of %d'%(PairN,len(Pairs)))
                FPrimer=PP.ForwardPrimer
                RPrimer=PP.ReversePrimer
                PP.FalseAmplicons=[]
                if FPrimer.SpecPenalty > SpecPenaltyThreshold and RPrimer.SpecPenalty <= SpecPenaltyThreshold:
                    PP.Penalty = PP.BasicPenalty + GenPar.SpecTestMaxPenalty
                    PP.ObjectiveSpecTest = False
                    PP.NonObjectiveSpecTestInfo='Forward primer is disastrous non-specific'
                    continue
                elif FPrimer.SpecPenalty <= SpecPenaltyThreshold and RPrimer.SpecPenalty > SpecPenaltyThreshold:
                    PP.Penalty = PP.BasicPenalty + GenPar.SpecTestMaxPenalty
                    PP.ObjectiveSpecTest = False
                    PP.NonObjectiveSpecTestInfo='Reverse primer is disastrous non-specific'
                    continue
                elif FPrimer.SpecPenalty > SpecPenaltyThreshold and RPrimer.SpecPenalty > SpecPenaltyThreshold:
                    PP.Penalty = PP.BasicPenalty + 2*GenPar.SpecTestMaxPenalty
                    PP.NonObjectiveSpecTestInfo='Both primers are disastrous non-specific'
                    PP.ObjectiveSpecTest = False
                    continue

                for ChrName in ChrSeqs.keys():
                    if (len(FPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName]) + len(RPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName])) * (len(FPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName]) + len(RPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName])) >= 500000000:
                        PP.Penalty = PP.BasicPenalty + GenPar.SpecTestMaxPenalty
                        if PP.ObjectiveSpecTest:
                            PP.NonObjectiveSpecTestInfo='Primers are disastrous non-specific at chr. %s'%ChrName
                        else:
                            PP.NonObjectiveSpecTestInfo=', %s'%ChrName
                        PP.ObjectiveSpecTest = False
                        continue

                    elif PP.ObjectiveSpecTest == False:
                        continue

                    CurrentFalPos_Forward={}
                    CurrentFalPos_Reverse={}

                    FPr_FSt = FPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName]
                    RPr_FSt = RPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName]
                    FPr_FSt_PosArray = FPrimer.FalsePositiveStarts_ForwardStrand_Sorted[ChrName]
                    RPr_FSt_PosArray = RPrimer.FalsePositiveStarts_ForwardStrand_Sorted[ChrName]

                    FPr_RSt = FPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName]
                    RPr_RSt = RPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName]
                    FPr_RSt_PosArray = FPrimer.FalsePositiveStarts_ReverseStrand_Sorted[ChrName]
                    RPr_RSt_PosArray = RPrimer.FalsePositiveStarts_ReverseStrand_Sorted[ChrName]

                    "SEARCH FOR POSSIBLE FALSE_POSITIVE AMPLICONS"

                    ExploringDist = ParDesired.MinFalseAmpliconSize
                    for FwdMMcounts,FwdPosArray,RevMMcounts,RevPosArray in ([FPr_FSt,FPr_FSt_PosArray,FPr_RSt,FPr_RSt_PosArray],
                            [RPr_FSt,RPr_FSt_PosArray,RPr_RSt,RPr_RSt_PosArray],
                            [FPr_FSt,FPr_FSt_PosArray,RPr_RSt,RPr_RSt_PosArray],
                            [RPr_FSt,RPr_FSt_PosArray,FPr_RSt,FPr_RSt_PosArray]):
                        if len(FwdPosArray) <= len(RevPosArray):      # выбор того массива в качестве рэперного, в котором меньше элементов
                            for FalPosFwd in FwdPosArray:
                                if len(RevPosArray)==0:
                                    continue
                                Range = GetRange(RevPosArray,FalPosFwd-ExploringDist-20,FalPosFwd+ExploringDist+20)
                                if Range.StartIndex == None:
                                    continue
                                for FalPosRev in RevPosArray[Range.StartIndex:Range.EndIndex+1]:
                                    if FalPosRev - FalPosFwd > 5 and FalPosRev - FalPosFwd < ParDesired.MinFalseAmpliconSize:
                                        FalseAmp=TFalseAmplicon()
                                        FalseAmp.Chr=ChrName
                                        FalseAmp.Start=FalPosFwd
                                        FalseAmp.End=FalPosRev
                                        FalseAmp.Size=FalPosRev-FalPosFwd+1
                                        FalseAmp.ForwardPrimerMismatches=FwdMMcounts[FalPosFwd]
                                        FalseAmp.ReversePrimerMismatches=RevMMcounts[FalPosRev]
                                        FalseAmp.Penalty = (((ParDesired.MinFalseAmpliconSize*1.15 - FalseAmp.Size)/100)**2.5)*(7 - FalseAmp.ForwardPrimerMismatches - FalseAmp.ReversePrimerMismatches )
                                        PP.FalseAmplicons.append(FalseAmp)
                        else:
                            for FalPosRev in RevPosArray:
                                Range = GetRange(FwdPosArray,FalPosRev-ExploringDist-20,FalPosRev+ExploringDist+20)
                                if Range.StartIndex == None:
                                    continue
                                for FalPosFwd in FwdPosArray[Range.StartIndex:Range.EndIndex+1]:
                                    if FalPosRev - FalPosFwd > 5 and FalPosRev - FalPosFwd < ParDesired.MinFalseAmpliconSize:
                                        FalseAmp=TFalseAmplicon()
                                        FalseAmp.Chr=ChrName
                                        FalseAmp.Start=FalPosFwd
                                        FalseAmp.End=FalPosRev
                                        FalseAmp.Size=FalPosRev-FalPosFwd+1
                                        FalseAmp.ForwardPrimerMismatches=FwdMMcounts[FalPosFwd]
                                        FalseAmp.ReversePrimerMismatches=RevMMcounts[FalPosRev]
                                        FalseAmp.Penalty = (((ParDesired.MinFalseAmpliconSize*1.15 - FalseAmp.Size)/100)**2.5)*(7 - FalseAmp.ForwardPrimerMismatches - FalseAmp.ReversePrimerMismatches )
                                        PP.FalseAmplicons.append(FalseAmp)

                    # CurrentFalPos_Forward={}
                    # CurrentFalPos_Forward.update(FPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName])
                    # CurrentFalPos_Forward.update(RPrimer.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand[ChrName])
                    #
                    # CurrentFalPos_Reverse={}
                    # CurrentFalPos_Reverse.update(FPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName])
                    # CurrentFalPos_Reverse.update(RPrimer.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand[ChrName])
                    #
                    # for FalPosFwd in CurrentFalPos_Forward.keys():
                    #     for FalPosRew in CurrentFalPos_Reverse.keys():
                    #         if FalPosRew - FalPosFwd > 5 and FalPosRew - FalPosFwd < SearchParameters.MaxFalsePCR_ProductSize:
                    #             FalseAmp=TFalseAmplicon()
                    #             FalseAmp.Chr=ChrName
                    #             FalseAmp.Start=FalPosFwd
                    #             FalseAmp.End=FalPosRew
                    #             FalseAmp.Size=FalPosRew-FalPosFwd+1
                    #             FalseAmp.ForwardPrimerMismatches=CurrentFalPos_Forward[FalPosFwd]
                    #             FalseAmp.ReversePrimerMismatches=CurrentFalPos_Reverse[FalPosRew]
                    #             FalseAmp.Penalty = (((5000 - FalseAmp.Size)/100)**2.5)*(7 - FalseAmp.ForwardPrimerMismatches - FalseAmp.ReversePrimerMismatches )
                    #             PP.FalseAmplicons.append(FalseAmp)
                    #

                PP.CalculateFullPenalty(ParDesired,ParPenalty,GenPar,True,Compl)

                PP.FalseAmplicons.sort(key=GetPenalty,reverse=True)




            if len(Pairs) < 10000:
                PairDistributionByFalseAmpliconsCount={}
                for PP in Pairs:
                    try:
                        PairDistributionByFalseAmpliconsCount[len(PP.FalseAmplicons)]+=1
                    except:
                        PairDistributionByFalseAmpliconsCount[len(PP.FalseAmplicons)]=1

                AuxInfoFile.write('\n\nDistribution of PAIRS by false amplicons count:')
                for K in sorted(PairDistributionByFalseAmpliconsCount.keys()):
                    AuxInfoFile.write('\n%d f.a.: %d pairs'%(K,PairDistributionByFalseAmpliconsCount[K]))


            "DELETING REDUNDANT AND NON_NEEDED INFO IN ORDER TO REDUCE DUMP SIZE ON DISK"

            for P in Pairs:
                P.CompleteFalseAmpliconsCount = len(P.FalseAmplicons)
                P.FalseAmplicons = P.FalseAmplicons[:GenPar.FalseAmpliconsToShow*2]

            for Pr in SpecTestForwardPrimers + SpecTestReversePrimers:
                del Pr.FalsePositiveStarts_ForwardStrand_Sorted
                del Pr.FalsePositiveStarts_ReverseStrand_Sorted
                del Pr.MismatchCountsByChrByFalsePositiveStarts_ForwardStrand
                del Pr.MismatchCountsByChrByFalsePositiveStarts_ReverseStrand
                del Pr.FalsePositivePositions

            print('\nDumping...')
            with open(PickleDB_FileNameBase+'.level2.db','wb') as output:
                pickle.dump((SpecPenaltyThreshold,SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID), output, pickle.HIGHEST_PROTOCOL)
            print('Completed...')




        if opts.OnlyVisualizationOfTCGAData:
            CompletePairsList = []
        elif DumpExists[2]:
            print('\nLoading from file %s...'%(PickleDB_FileNameBase+'.level3.db'))
            with open(PickleDB_FileNameBase+'.level3.db','rb',10000000) as input:
                CompletePairsList,SpecPenaltyThreshold,SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID = pickle.load(input)
            print('Completed.')
        else:
            # print('test1')
            #######################################№№№№№№№№№######################################################################################
            # Составление композиции нескольких пар - МЕТОД 1 (СТАРЫЙ)      ##############################################################################
            ################################################№№####################################################################################

            CompositionMethod='NewMultiThread'

            if CompositionMethod=='Old':
                class TUncoveredArea():
                    def __init__(self):
                        self.Start=0
                        self.End=0

                MajorIterationNumber = 0
                MaxMajorIterationNumber = 4

                UncoveredAreas=[]
                InitUA=TUncoveredArea()
                InitUA.Start=start+SideGap
                InitUA.End=fin-SideGap
                UncoveredAreas.append(InitUA)

                PairsSequencesByIteration=[]
                for MajorIterationNumber in range(MaxMajorIterationNumber):
                    print('\n\n###############################################\n###############################################')
                    AreaDescriptions=[]
                    for UA in UncoveredAreas:
                        AreaDescriptions.append('%d [%d-%d]'%(UA.End-UA.Start+1,UA.Start,UA.End))
                    print('Iteration %d. Area sizes to be covered by amplicons: %s\n'%(MajorIterationNumber+1, '; '.join( AreaDescriptions) ))
                    UA_ForThisIteration=[]
                    NewUncoveredAreas=[]   # Непокрытые ампликонами участки, которые появляются на НОВОЙ итерации
                    PairsSequence=[]
                    # CurrentIterationUA_Starts=sorted(UncoveredAreasByUA_Starts.keys())
                    for UA in UncoveredAreas:
                        NewUA_Start = None
                        print('\n')
                        print('>>>Processing area of %d bp [%d-%d]'%(UA.End-UA.Start+1,UA.Start,UA.End))
                        # UncoveredAreasByUA_Starts[St]

                        # if MajorIterationNumber==1:
                        #     Pause=1

                        CurrentUA_Start=UA.Start
                        CurrentUA_End=UA.End

                        PenaltyIncreaseStep=23
                        PenaltyThreshold=12
                        InitialPenaltyThreshold=12
                        MaxDrawbackToIncreasePenaltyInitial=ParLimit.MaxAmpliconLength/7+50 # ИСХОДНОЕ разрешённое over-перекрытие двух ампликонов, чтобы ещё не увеличивать порог штрафов
                        MaxAllowedPenaltyBeforeIncreaseMaxDrawbackLimit = 170+(MajorIterationNumber*40)        # максимально достигнутый порог по штрафам, прежде чем увеличиать Drawback
                        MaxDrawbackToIncreasePenaltyCurrent=ParLimit.MaxAmpliconLength/5+50  # ТЕКУЩЕ разрешённое перекрытия двух ампликонов, чтобы ещё не увеличивать порог штрафов
                        OverDrawbackIncreaseStep=int(MaxDrawbackToIncreasePenaltyInitial/4)
                        # MaxAllowedOverDrawback=MaxDrawbackToIncreasePenaltyInitial*2      # максимально вообще разрешённый порог Drawbacka
                        MaxAllowedOverDrawback=min(MaxDrawbackToIncreasePenaltyInitial*(2.2+0.15*MajorIterationNumber),ParLimit.MinAmpliconLength/1.4)      # максимальнр вообще разрешённый порог Drawbacka

                        EmergenceStepover=25  # когда ничего не помогает, тогда сдвигаемся направо на 100 п.н.

                        CurrentScanEnd=max(CurrentUA_Start - int(ParLimit.MaxPrimerLength),30)
                        CurrentScanStart=max(0,CurrentScanEnd-MaxDrawbackToIncreasePenaltyCurrent)

                        AtLeastOnePairAssigned=False
                        while True:
                            OptimalPair=None
                            OptimalUpCoverage=15
                            for P in Pairs:
                                if P.Penalty > PenaltyThreshold:
                                    continue
                                if P.ForwardPrimer.InChrStart <= CurrentScanEnd and P.ForwardPrimer.InChrStart >= CurrentScanStart:
                                    if len(PairsSequence):
                                        CurrentUpCoverage = min(UA.End, P.ReversePrimer.InChrEnd) - max(CurrentScanEnd,PairsSequence[-1].ReversePrimer.End-int(MinAmpliconOverlap/2))
                                    else:
                                        CurrentUpCoverage = min(UA.End, P.ReversePrimer.InChrEnd) - max(CurrentScanEnd,UA.Start)
                                    # if
                                    if OptimalUpCoverage < CurrentUpCoverage:
                                        OptimalUpCoverage=CurrentUpCoverage
                                        OptimalPair=P
                            if OptimalPair!=None:
                                print('Pair %d-%d with penalty %.1f has been added to the sequence (%d upcoverage)'%(OptimalPair.ForwardPrimer.InChrStart,OptimalPair.ReversePrimer.InChrEnd,OptimalPair.Penalty,OptimalUpCoverage))
                                AtLeastOnePairAssigned=True
                                if NewUA_Start!=None:
                                    NewUA_End = OptimalPair.ForwardPrimer.InChrStart
                                    NewUA=TUncoveredArea()
                                    NewUA.Start = NewUA_Start
                                    NewUA.End = NewUA_End
                                    if NewUA.End > NewUA.Start+5:
                                        NewUncoveredAreas.append(NewUA)
                                    NewUA_Start=None

                                PairsSequence.append(OptimalPair)
                                CurrentScanEnd=OptimalPair.ReversePrimer.InChrEnd-MinAmpliconOverlap
                                CurrentScanStart=max(CurrentScanStart+ParLimit.MaxPrimerLength+5,max(0,CurrentScanEnd-MaxDrawbackToIncreasePenaltyCurrent))
                            else:
                                PenaltyThreshold+=PenaltyIncreaseStep
                                if PenaltyThreshold > MaxAllowedPenaltyBeforeIncreaseMaxDrawbackLimit:

                                    PenaltyThreshold = InitialPenaltyThreshold
                                    MaxDrawbackToIncreasePenaltyCurrent += OverDrawbackIncreaseStep
                                    if  MaxDrawbackToIncreasePenaltyCurrent > MaxAllowedOverDrawback:
                                        print('appropriate pair with start located at %d-%d cannot be found. moving forward to %d bp'%(CurrentScanStart,CurrentScanEnd,EmergenceStepover))
                                        if MajorIterationNumber==1:
                                            Pause=0
                                        if NewUA_Start==None:
                                            if len(PairsSequence)>0:
                                                NewUA_Start = PairsSequence[-1].ReversePrimer.InChrEnd
                                            else:
                                                NewUA_Start = UA.Start
                                        CurrentScanEnd += EmergenceStepover  # увеличение диапазона поиска пар (передний, правй конец). когда ничего не помогает
                                        MaxDrawbackToIncreasePenaltyCurrent=MaxDrawbackToIncreasePenaltyInitial
                                    else:
                                        print('Increasing drawback to %d'%MaxDrawbackToIncreasePenaltyCurrent)

                                    CurrentScanStart=max(0,CurrentScanEnd-MaxDrawbackToIncreasePenaltyCurrent)  # увеличение диапазона поиска пар (левый конец)

                                # MaxDrawbackToIncreasePenaltyCurrent=MaxDrawbackToIncreasePenaltyInitial
                            if CurrentScanStart > CurrentUA_End-SideGap+int(ParLimit.MaxPrimerLength+4):
                                print('Total %d pairs (across this iteration)'%len(PairsSequence))
                                break

                        if not AtLeastOnePairAssigned:
                            NewUA=TUncoveredArea()
                            NewUA.Start = UA.Start
                            NewUA.End = UA.End
                            if NewUA.End > NewUA.Start+5:
                                NewUncoveredAreas.append(NewUA)

                        else:
                            if NewUA_Start!=None:
                                NewUA_End = UA.End
                                NewUA=TUncoveredArea()

                                NewUA.Start = NewUA_Start
                                NewUA.End = NewUA_End
                                if NewUA.End > NewUA.Start+5:
                                    NewUncoveredAreas.append(NewUA)

                    UncoveredAreas=NewUncoveredAreas

                    PairsSequencesByIteration.append(PairsSequence)

                    if len(UncoveredAreas)==0:
                        print('Sequence coverage COMPLETED.\nTotal iterations: %d. Total pairs by iterations: %s'%(MajorIterationNumber+1, str(   [len(x) for x in PairsSequencesByIteration   ] ) ))
                        break


                if len(UncoveredAreas)>0:
                    print('Sequence coverage INCOMPLETED.\nTotal iterations: %d. Total pairs by iterations: %s'%(MajorIterationNumber+1, str(   [len(x) for x in PairsSequencesByIteration   ] ) ))
                    print('Remaining areas: %s'%('; '.join( sorted([str(x.End-x.Start+1) for x in UncoveredAreas ]) ) ))
                    UncoveredAreasSequences = [ ChrSeqs[Chr].seq[x.Start:x.End+1] for x in UncoveredAreas ]
                    SumUncoveredLength = sum( [len(x) for x in  UncoveredAreasSequences ]  )
                    SumUncoveredCpG=sum( [x.count('CG') for x in  UncoveredAreasSequences ]  )
                    print('Total not-covered sequence length: %d with %d CpG dinucleotides' %(SumUncoveredLength,SumUncoveredCpG))

                CompletePairsList=[]
                for PS in PairsSequencesByIteration:
                    for P in PS:
                        CompletePairsList.append(P)

            #######################################№№№№№№№№№######################################################################################
            # Составление композиции нескольких пар - МЕТОД 2 (НОВЫЙ + МНОГОПОТОЧНЫЙ)      ########################################################
            ################################################№№####################################################################################
            elif CompositionMethod=='NewMultiThread':
                CompletePairsList=[]
                ThreadsCount=GenPar.CPUcount
                # GenPar.AttemptsCountPerThread   # число попыток, в каждой из которых происходит случайный выбор
                                         # пары и дальнейшее наращивание последовательности ампликонов вправо и влево

                # if
                SourceCpG_coords=[]
                Pos=0
                while True:         # Составление списка координат CpG-динуклеотидов
                    Pos=SourceSeq.find('CG',Pos)
                    if Pos > 0:
                        SourceCpG_coords.append(Pos)
                    else:
                        break
                    Pos+=1


                PairsMT=[]
                PairsByShortID={}
                CurrentPair_ID=0
                for P in Pairs:
                    P.ShortID = CurrentPair_ID
                    PairsByShortID[CurrentPair_ID]=P
                    PairMT = TPrimerPair()
                    PairMT.CopySkeletonFrom(P)
                    PairsMT.append(PairMT)
                    CurrentPair_ID+=1

                pool = Pool(processes=ThreadsCount)
                StartupArray=[]
                for n in range(ThreadsCount):
                    ComPar = TCompositePairsStartupParameters()
                    ComPar.SourceSeq=SourceSeq
                    ComPar.AttemptsPerThread=int(GenPar.TotalBootstraps/GenPar.CPUcount+0.5)
                    ComPar.Pairs=PairsMT
                    ComPar.SourceCpG_coords=SourceCpG_coords
                    ComPar.SideGap=SideGap
                    ComPar.ThreadNumber=n
                    ComPar.Mode=opts.Mode
                    StartupArray.append(ComPar)
                ResultsArray = pool.map(CompositePairs, StartupArray)
                pool.close()
                pool.join()


                CompletePairsID_ListByThread=[]
                BestScoresByThread=[]
                for ThreadN in range(ThreadsCount):
                    CurrentCompletePairsID_List,CurrentBestScore = ResultsArray[ThreadN]
                    CompletePairsID_ListByThread.append(CurrentCompletePairsID_List)
                    BestScoresByThread.append(CurrentBestScore)

                BestScore=None
                for N in range(len(CompletePairsID_ListByThread)):
                    if BestScore == None:
                        BestScore = BestScoresByThread[N]
                        CompletePairsID_List = CompletePairsID_ListByThread[N]
                    if BestScore < BestScoresByThread[N]:
                        BestScore = BestScoresByThread[N]
                        CompletePairsID_List = CompletePairsID_ListByThread[N]

                CompletePairsList = [PairsByShortID[x] for x in CompletePairsID_List]



                # BestAttemptN = ScoresByAttempt.index(max(ScoresByAttempt))
                # CompletePairsList = PairSequencesByAttempt[BestAttemptN]


            print('Dumping...')
            with open(PickleDB_FileNameBase+'.level3.db','wb') as output:
                pickle.dump((CompletePairsList,SpecPenaltyThreshold,SpecTestForwardPrimers,SpecTestReversePrimers,ForwardPrimers,ReversePrimers,DesiredForwardPrimers,DesiredReversePrimers,BulkForwardPrimers,BulkReversePrimers,PrimersByID,Pairs,PairsByID), output, pickle.HIGHEST_PROTOCOL)
            print('Completed...')



































    #   Рисование расположения праймеровъ и всего прочего


        ImageFileName=None
        # if 'win' in sys.platform:

        MainCanvasW=1400
        MainCanvasH=200
        SupplCanvasW=MainCanvasW
        SupplCanvasH=75
        ENCODE_state_Canvas_H = SupplCanvasH/2
        CG_profileH=100
        CoverageProfileH=100
        MapW=MainCanvasW
        MapH=int(MainCanvasH*0.5)
        MapWStart=int((MainCanvasW-MapW)/2)
        MapHStart=int((MainCanvasH-MapH)/2)
        LargeMap=[0]*10*len(SourceSeq)
        LargeMapSize=10*len(SourceSeq)

        Pos=0
        OverDraw = max(1, int(1+ (fin - start - 2000)/3500 + 0.5 ))
        CpG_Count=0
        while Pos>=0:
            Pos=SourceSeq.find('CG',Pos+1)
            if Pos>0:
                CpG_Count+=1
                LargeMapStart=max(10*Pos-OverDraw,0)
                LargeMapEnd=min(10*(Pos+1)+OverDraw,LargeMapSize)
                for i in range(LargeMapStart,LargeMapEnd):
                    LargeMap[i]+=1

        SmallMap=[0]*MapW
        for i in range(MapW):
            LargeMapStart=int(i/MapW * LargeMapSize)
            LargeMapEnd=int((i+1)/MapW * LargeMapSize)
            SmallMap[i]=sum(LargeMap[LargeMapStart:LargeMapEnd]) / (LargeMapEnd-LargeMapStart+0.1)

#  загрузка данныъ TCGA из кэша на диске
#         TCGA_db_FileName = os.path.join(opts.PrimerDatabaseFolder,'%s.%s.%s.%s.TCGA_info.db'%(Gene,ParametersHash,SeqHash,hashmd5(opts.CrossHubTCGAMethylationResultsFiles)))
#         if os.path.exists(TCGA_db_FileName) and not opts.UpdateTCGA_Cache:
#             with open(TCGA_db_FileName,'rb') as input:
#                 CurrentCpGsAcrossAllFiles,TCGA_ResultsFileDesignations = pickle.load(input)
#             print('\nCpG methylation TCGA for gene %s info is loaded from cache'%(Gene))
#         else:
#             CurrentCpGsAcrossAllFiles = {}
#             for Pos in range(start,fin+1):
#                 try:
#                     CpG = CpGsByChrByPos[Chr][Pos]
#                 except KeyError:
#                     continue
#                 try:
#                     CurrentCpGsAcrossAllFiles[Chr]
#                 except KeyError:
#                     CurrentCpGsAcrossAllFiles[Chr] = {}
#                 CurrentCpGsAcrossAllFiles[Chr][Pos] = CpG
#             with open(PickleDB_FileNameBase+'.TCGA_info.db','wb') as output:
#                 pickle.dump((CurrentCpGsAcrossAllFiles,TCGA_ResultsFileDesignations), output, pickle.HIGHEST_PROTOCOL)
#                 print('\nCpG methylation TCGA for gene %s info is cached'%(Gene))

#  определяем, существуют ли данные о коррелциях между метилированием и РНАсек для разных локализаций
#         CorrelationContainingTCGA_FilesCount = 0
#         for FileDesignation in TCGA_ResultsFileDesignations:
#             CurrentCpGs = []
#             for Pos in range(start,fin+1):
#                 try:
#                     CpG = CurrentCpGsAcrossAllFiles[Chr][Pos]
#                 except KeyError:
#                     continue
#                 try:
#                     CpG.Rs_paired[FileDesignation]
#                     CpG.Rs_pooled[FileDesignation]
#                 except KeyError:
#                     continue
#                 CurrentCpGs.append(CpG)
#             if len(CurrentCpGs) != 0:
#                 CorrelationContainingTCGA_FilesCount += 1


#  размечаем координаты различных объектов

        # Canvas=Canvas(bg="white", width=MainCanvasW, height=MainCanvasH+SupplCanvasH*len(CompletePairsList))

        if opts.OnlyVisualizationOfTCGAData:
            ENCODE_SegDescriptions_OffsetH = MainCanvasH
        else:
            ENCODE_SegDescriptions_OffsetH = MainCanvasH  +CG_profileH*8+SupplCanvasH*len(CompletePairsList)

        TCGA_ResultsOffset_T = ENCODE_SegDescriptions_OffsetH + ENCODE_state_Canvas_H*len(SegDescriptions) + 40
        TCGA_ColorSchemaOffsetH_T = TCGA_ResultsOffset_T   + SupplCanvasH*len(CrossHub_data_by_TCGA_code) + 20
        if opts.MergeTCGAHyperAndHypoMethylation:
            SecondMainCanvasOffsetH = TCGA_ColorSchemaOffsetH_T + SupplCanvasH
        else:
            TCGA_ResultsOffset_C      = TCGA_ColorSchemaOffsetH_T + SupplCanvasH + 20
            TCGA_ColorSchemaOffsetH_C = TCGA_ResultsOffset_C   + SupplCanvasH*len(CrossHub_data_by_TCGA_code) + 20
            SecondMainCanvasOffsetH = TCGA_ColorSchemaOffsetH_C + SupplCanvasH
        Complete_ImageH = int(SecondMainCanvasOffsetH + MainCanvasH - 30)

        IImage=Image.new('RGB', (MainCanvasW,Complete_ImageH),(255,255,255))
        IDraw=ImageDraw.Draw(IImage)

## Drawing main map with CpG sites location
        MaxS = 1.5
        for i in range(MapW):

            SCount = min(SmallMap[i],MaxS)

            RColorBase1=255
            GColorBase1=238
            BColorBase1=169
            RColorBase2=0
            GColorBase2=126
            BColorBase2=255
            RColor=int(RColorBase1+(RColorBase2-RColorBase1)*SCount/MaxS)
            GColor=int(GColorBase1+(GColorBase2-GColorBase1)*SCount/MaxS)
            BColor=int(BColorBase1+(BColorBase2-BColorBase1)*SCount/MaxS)
            CColor=color(RColor,GColor,BColor)

            # Canvas.create_line(coord,fill=CColor)
            coord=i,MapHStart,i,MapHStart+MapH
            IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)
            coord= i,MapHStart+SecondMainCanvasOffsetH,i,MapHStart+MapH*0.2+SecondMainCanvasOffsetH
            IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)

#  отображение шкалы координат
        IFont = ImageFont.truetype("verdana.ttf", size=13)
        TickedCurrentTickGenomeCoords = []
        # if fin - start
        for TickInterval,TickSize in zip([100000,50000,10000,5000,1000,500,100],[0.21,0.21,0.17,0.14,0.11,0.07,0.04]):
            TickNumber = 0
            # print([TickInterval,TickSize])
            ImageTickInterval = MapW * TickInterval / (fin - start +1)
            if ImageTickInterval < 10:
                continue
            while True:
                CurrentTickGenomeCoord = start + (TickInterval - (start%TickInterval)) + TickInterval*TickNumber
                if CurrentTickGenomeCoord > fin:
                    break
                if CurrentTickGenomeCoord in TickedCurrentTickGenomeCoords:
                    TickNumber+=1
                    continue
                TickedCurrentTickGenomeCoords.append(CurrentTickGenomeCoord)
                CurrentTickImageCoord = (CurrentTickGenomeCoord - start)/(fin - start) * MapW
                coord = CurrentTickImageCoord,MapHStart+MapH*0.22+SecondMainCanvasOffsetH,CurrentTickImageCoord,MapHStart+MapH*(0.22+TickSize)+SecondMainCanvasOffsetH
                IDraw.line(coord,fill=(30,30,30),width=1)
                if TickInterval > (fin - start)/4:
                    TickText = intWithCommas(CurrentTickGenomeCoord)
                else:
                    if CurrentTickGenomeCoord%1000 > 0:
                        TickText = intWithCommas(CurrentTickGenomeCoord%1000)
                    else:
                        TickText = intWithCommas(CurrentTickGenomeCoord)

                w,h = IDraw.textsize(text=TickText,font=IFont)
                if w * (fin - start)/TickInterval < 0.7 * MapW:
                    IDraw.text(xy=(max(CurrentTickImageCoord-w/2,10), MapHStart+MapH*(0.22+TickSize + 0.01)+SecondMainCanvasOffsetH),text=TickText,fill=(0,0,0),font=IFont)
                TickNumber+=1


#       рисование подписи, что это за ген и где начинается (вверху и внизу)
        # Canvas.create_text(20, 10, anchor=NW,font="Purisa")
        IDraw.text(xy=(20,10),text=intWithCommas(start),fill=(0,0,0),font=IFont)
        IDraw.text(xy=(20,10+SecondMainCanvasOffsetH),text=intWithCommas(start),fill=(0,0,0),font=IFont)


        # Canvas.create_text(20, 27, anchor=NW,font="Purisa",text=Gene)
        # w,h = IDraw.textsize(text=Gene,font=IFont)
        IDraw.text(xy=(20,27),text=Gene,fill=(0,0,0),font=IFont)
        IDraw.text(xy=(20,27+SecondMainCanvasOffsetH),text=Gene,fill=(0,0,0),font=IFont)

        # Canvas.create_text(MainCanvasW-20, 10, anchor=NE,font="Purisa",text=intWithCommas(fin))
        w,h = IDraw.textsize(text=intWithCommas(fin),font=IFont)
        IDraw.text(xy=(MainCanvasW-20-w, 10),text=intWithCommas(fin),fill=(0,0,0),font=IFont)
        IDraw.text(xy=(MainCanvasW-20-w, 10+SecondMainCanvasOffsetH),text=intWithCommas(fin),fill=(0,0,0),font=IFont)


        text='%s bp., %d CpG'%(intWithCommas(fin-start),CpG_Count)
        # Canvas.create_text(MainCanvasW-20, 27, anchor=NE,font="Purisa",text=text)
        w,h = IDraw.textsize(text=text,font=IFont)
        IDraw.text(xy=(MainCanvasW-20-w, 27),text=text,fill=(0,0,0),font=IFont)
        IDraw.text(xy=(MainCanvasW-20-w, 27+SecondMainCanvasOffsetH),text=text,fill=(0,0,0),font=IFont)

        coord = int(MapWStart + MapW*SideGap/len(SourceSeq) ),MapHStart+MapH+10, int(MapWStart+MapW - MapW*SideGap/len(SourceSeq) ),MapHStart+MapH+15
        # Canvas.create_rectangle(coord,fill='#c1f1a7',width=0)
        IDraw.rectangle(xy=coord,fill=(193,241,167))

        coord = int(MapWStart + MapW*SideGap/len(SourceSeq) ),MapHStart+MapH*0.5+10+SecondMainCanvasOffsetH, int(MapWStart+MapW - MapW*SideGap/len(SourceSeq) ),MapHStart+MapH*0.5+15+SecondMainCanvasOffsetH
        IDraw.rectangle(xy=coord,fill=(193,241,167))

########################################################################
##  ENCODE chromatin state segmentation info
########################################################################
        if len(SegDescriptions) > 0:
            for n in range(len(SegDescriptions)):
                SegDes = SegDescriptions[n]
                for i in range(MapW):
                    genome_coord = int(start + float(fin - start)*float(i)/float(MapW))
                    # print(genome_coord)
                    SegType = RevealSegTypeByChrByPosition(SegDes.SegListsByChr,Chr,genome_coord)
                    if SegType is None: continue
                    if SegType in ColorsBy_ENCODE_SegState:
                        CColor = ColorsBy_ENCODE_SegState[SegType]
                    else: CColor = '#DDDDDD'

                    coord=i,ENCODE_SegDescriptions_OffsetH + ENCODE_state_Canvas_H*(n + 0.45),i,ENCODE_SegDescriptions_OffsetH + ENCODE_state_Canvas_H*(n + 0.9)
                    IDraw.line(coord,fill=CColor,width=1)
                ENCODE_text_descr = "ENCODE %s %s"%(SegDes.CellType,SegDes.ML_Type)
                IDraw.text(xy=(20, ENCODE_SegDescriptions_OffsetH + ENCODE_state_Canvas_H*(n)),text=ENCODE_text_descr,fill=(0,0,0),font=IFont)


########################################################################
#  Visualizing TCGA Methylation
####################################################################

        if len(CrossHub_data_by_TCGA_code) > 0:
    #       рисование плохих и хороших CpG-положений
            if opts.MergeTCGAHyperAndHypoMethylation:
    #       общая цветовая схема - HYPERMETHYLATION
                coord = MapW*0.2+20, TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.25,MapW*0.2+35,TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.45
                IDraw.rectangle(xy=coord,fill=(225,225,225))
                for i in range (int(MapW*0.2)):
                    GradientCoord = i/(MapW*0.2)
                    GradientCoord = GradientCoord**0.74
                    RGBColor = Gradient(HypermethylationGradientScorePoints,GradientCoord)
                    coord = 35 + MapW*0.2 + i, TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.25,35 + MapW*0.2 + i,TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.45
                    IDraw.rectangle(xy=coord,fill=RGBColor)

                    RGBColor = Gradient(HypomethylationGradientScorePoints,GradientCoord)
                    coord = 20 + MapW*0.2 - i, TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.25,20 + MapW*0.2 - i,TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.45
                    IDraw.rectangle(xy=coord,fill=RGBColor)

                IDraw.text(xy=(17, TCGA_ColorSchemaOffsetH_T),text='100     Hypo-M score',fill=(0,0,0),font=IFont)
                IDraw.text(xy=(17 + MapW*0.2, TCGA_ColorSchemaOffsetH_T),text=' 0',fill=(0,0,0),font=IFont)
                IDraw.text(xy=(MapW*0.4 - 100, TCGA_ColorSchemaOffsetH_T),text=' Hyper-M score    100',fill=(0,0,0),font=IFont)
            else:
    #       общая цветовая схема - HYPERMETHYLATION
                coord = 20, TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.25,35,TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.45
                IDraw.rectangle(xy=coord,fill=(225,225,225))
                for i in range (int(MapW*0.25)):
                    GradientCoord = i/(MapW*0.25)
                    GradientCoord = GradientCoord**0.74
                    RGBColor = Gradient(HypermethylationGradientScorePoints,GradientCoord)

                    coord = 35 + i, TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.25,35+i,TCGA_ColorSchemaOffsetH_T + SupplCanvasH*0.45
                    IDraw.rectangle(xy=coord,fill=RGBColor)

                IDraw.text(xy=(17, TCGA_ColorSchemaOffsetH_T),text='0',fill=(0,0,0),font=IFont)
                IDraw.text(xy=(17 + MapW*0.25, TCGA_ColorSchemaOffsetH_T),text='100     TCGA score (hyper-meth)',fill=(0,0,0),font=IFont)


    #       общая цветовая схема - HYPoMETHYLATION
                coord = 20, TCGA_ColorSchemaOffsetH_C + SupplCanvasH*0.25,35,TCGA_ColorSchemaOffsetH_C + SupplCanvasH*0.45
                IDraw.rectangle(xy=coord,fill=(225,225,225))
                for i in range (int(MapW*0.25)):
                    GradientCoord = i/(MapW*0.25)
                    GradientCoord = GradientCoord**0.74
                    RGBColor = Gradient(HypermethylationGradientScorePoints,GradientCoord)
                    coord = 35 + i, TCGA_ColorSchemaOffsetH_C + SupplCanvasH*0.25,35+i,TCGA_ColorSchemaOffsetH_C + SupplCanvasH*0.45
                    IDraw.rectangle(xy=coord,fill=RGBColor)

                IDraw.text(xy=(17, TCGA_ColorSchemaOffsetH_C),text='0',fill=(0,0,0),font=IFont)
                IDraw.text(xy=(17 + MapW*0.25, TCGA_ColorSchemaOffsetH_C),text='100     TCGA score (hypo-meth)',fill=(0,0,0),font=IFont)


#       непосредственно рисунок плохих и хороших CpG ПО ОТДЕЛЬНОСТИ ДЛЯ КАЖДОГО ФАЙЛА -  HYPERMETHYLATION

        CpG_SmallMapsOverdraw = 2
        TCGA_DrawN = 0
        for TCGA_code in CrossHub_TCGA_codes:
            CrossHub_data = CrossHub_data_by_TCGA_code[TCGA_code]
            CurrentCpGs = []
            CurrentCpG_Positions = []
        # вытаскиваем CpG для интересующего участка
            CorrelationsPresent = False
            for Pos in range(start,fin+1):
                if not Chr in CrossHub_data.CpG_BasicInfo_ByChr_ByPos: continue
                if not Pos in CrossHub_data.CpG_BasicInfo_ByChr_ByPos[Chr]: continue

                CpG = CrossHub_data.CpG_BasicInfo_ByChr_ByPos[Chr][Pos]

                if len(CpG.Rs_paired) > 0 or len(CpG.Rs_pooled) > 0:
                    CorrelationsPresent = True

                CurrentCpGs.append(CpG)
                CurrentCpG_Positions.append(Pos)

            if opts.MergeTCGAHyperAndHypoMethylation:
                CurrentCpGs.sort(key = (lambda x: max(x.HyperMeth_score,x.HypoMeth_score)))
            else:
                CurrentCpGs.sort(key = (lambda x: x.HyperMeth_score))

        # непосредственно рисуем картинку
            if opts.MergeTCGAHyperAndHypoMethylation:
                for CpG_n in range(len(CurrentCpGs)):
                    CpG = CurrentCpGs[CpG_n]
                    HyperScore = CpG.HyperMeth_score
                    HypoScore = CpG.HypoMeth_score
                    if math.isnan(HyperScore) and math.isnan(HypoScore): continue
                    elif not math.isnan(HyperScore) and math.isnan(HypoScore): Score = HyperScore
                    elif math.isnan(HyperScore) and not math.isnan(HypoScore): Score = HypoScore
                    else:  Score = max(HyperScore,HypoScore)
                    if Score < 0.1:
                        RGBColor=(225,225,225)
                    else:
                        GradientCoord = max(0,min(1,(Score - 0.1)/0.9))
                        GradientCoord = GradientCoord**0.74
                        if HyperScore >= HypoScore:
                            RGBColor = Gradient(HypermethylationGradientScorePoints,GradientCoord)
                        else:
                            RGBColor = Gradient(HypomethylationGradientScorePoints,GradientCoord)

                    i = (CurrentCpG_Positions[CpG_n] - start)/(fin - start + 1)*MapW
                    coord= i,TCGA_ResultsOffset_T + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH,\
                             i + CpG_SmallMapsOverdraw,\
                             TCGA_ResultsOffset_T + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH + SupplCanvasH*0.3
                    IDraw.rectangle(xy = coord,fill=RGBColor)
            else:
                for CpG_n in range(len(CurrentCpGs)):
                    CpG = CurrentCpGs[CpG_n]
                    if math.isnan(CpG.HyperMeth_score): continue
                    Score = CpG.HyperMeth_score
                    if Score < 0.1:
                        RGBColor=(225,225,225)
                    else:
                        GradientCoord = max(0,min(1,(Score - 0.1)/0.9))
                        GradientCoord = GradientCoord**0.74
                        RGBColor = Gradient(HypermethylationGradientScorePoints,GradientCoord)

                    i = (CurrentCpG_Positions[CpG_n] - start)/(fin - start + 1)*MapW
                    coord= i,TCGA_ResultsOffset_T + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH,\
                             i + CpG_SmallMapsOverdraw,\
                             TCGA_ResultsOffset_T + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH + SupplCanvasH*0.3
                    IDraw.rectangle(xy = coord,fill=RGBColor)

            TextOut = '%s:'%(CrossHub_data.SamplingInfo.Dataset_Abbreviation)
            if len(CurrentCpGs) > 0:
                TextOut += ' %d T, %d N samples (%d paired)'%(
                    CrossHub_data.SamplingInfo.Methyl_samples_count_T,
                    CrossHub_data.SamplingInfo.Methyl_samples_count_C,
                    CrossHub_data.SamplingInfo.Methyl_samples_count_paired
                )
            if Gene in CrossHub_data.RNASeq_BasicInfo_ByGeneName:
                if not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled) and not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired):
                    TextOut += ', RNA-Seq LogFC = %.1f / %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled,CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired)
                elif not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled):
                    TextOut += ', RNA-Seq LogFC pooled = %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled)
                elif not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired):
                    TextOut += ', RNA-Seq LogFC paired = %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired)
            TextOut += ' [%s]'%(CrossHub_data.SamplingInfo.Dataset_Name)
            IDraw.text(xy=(20,TCGA_ResultsOffset_T + 0.01*SupplCanvasH + TCGA_DrawN*SupplCanvasH),text=TextOut,fill=(0,0,0),font=IFont)
            # TCGA_DrawN += 1

            # рисование графика с корреляциями
            if CorrelationsPresent:
                Cor_plot_rel_size = 0.35
                CurrentH_Offset = TCGA_ResultsOffset_T + 0.6*SupplCanvasH + TCGA_DrawN*SupplCanvasH
                coord = 0,CurrentH_Offset,MapW,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size
                IDraw.rectangle(xy = coord,fill=(242,250,255))
                coord = 0,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,MapW,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2
                IDraw.rectangle(xy = coord,fill=(236,236,198))
                for CpG_n in range(len(CurrentCpGs)):
                    CpG = CurrentCpGs[CpG_n]
                    add_text = None
                    if Gene in CpG.Rs_paired:
                        Rs_paired = CpG.Rs_paired[Gene]
                    else:
                        Rs_paired = CpG.Rs_paired[sorted(CpG.Rs_paired.keys())[0]]
                        add_text = 'correlations are shown for %s gene'%(sorted(CpG.Rs_paired.keys())[0])

                    # Rs_paired = (fabs(Rs_paired)**0.75)*np.sign(Rs_paired)
                    if Gene in CpG.Rs_pooled:
                        Rs_pooled = CpG.Rs_pooled[Gene]
                    else:
                        Rs_pooled = CpG.Rs_pooled[sorted(CpG.Rs_pooled.keys())[0]]
                        add_text = 'correlations are shown for %s gene'%(sorted(CpG.Rs_pooled.keys())[0])
                    # Rs_pooled = (fabs(Rs_pooled)**0.75)*np.sign(Rs_pooled)

                    i = (CurrentCpG_Positions[CpG_n] - start)/(fin - start + 1)*MapW
                    if not math.isnan(Rs_pooled):
                        RGBColor=(68,120,251)
                        coord= i ,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,i,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2 - SupplCanvasH*Cor_plot_rel_size/2*Rs_pooled
                        IDraw.rectangle(xy = coord,fill=RGBColor)

                    if not math.isnan(Rs_paired):
                        RGBColor=(241,126,56)
                        coord= i+1 ,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,i+1,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2 - SupplCanvasH*Cor_plot_rel_size/2*Rs_paired
                        IDraw.rectangle(xy = coord,fill=RGBColor)
            TCGA_DrawN += 1




#       непосредственно рисунок плохих и хороших CpG ПО ОТДЕЛЬНОСТИ ДЛЯ КАЖДОГО ФАЙЛА -  HYPOMETHYLATION

        if not opts.MergeTCGAHyperAndHypoMethylation:
            CpG_SmallMapsOverdraw = 2
            TCGA_DrawN = 0
            for TCGA_code in CrossHub_TCGA_codes:
                CrossHub_data = CrossHub_data_by_TCGA_code[TCGA_code]
                CurrentCpGs = []
                CurrentCpG_Positions = []
                CorrelationsPresent = False
                for Pos in range(start,fin+1):
                    if not Chr in CrossHub_data.CpG_BasicInfo_ByChr_ByPos: continue
                    if not Pos in CrossHub_data.CpG_BasicInfo_ByChr_ByPos[Chr]: continue

                    CpG = CrossHub_data.CpG_BasicInfo_ByChr_ByPos[Chr][Pos]

                    if len(CpG.Rs_paired) > 0 or len(CpG.Rs_pooled) > 0:
                        CorrelationsPresent = True

                    CurrentCpGs.append(CpG)
                    CurrentCpG_Positions.append(Pos)

                CurrentCpGs.sort(key = (lambda x: x.HypoMeth_score))

                for CpG_n in range(len(CurrentCpGs)):
                    CpG = CurrentCpGs[CpG_n]
                    if math.isnan(CpG.HypoMeth_score): continue
                    Score = CpG.HypoMeth_score
                    if Score < 0.1:
                        RGBColor=(225,225,225)
                    else:
                        GradientCoord = max(0,min(1,(Score - 0.1)/0.9))
                        GradientCoord = GradientCoord**0.74
                        RGBColor = Gradient(HypomethylationGradientScorePoints,GradientCoord)

                    i = (CurrentCpG_Positions[CpG_n] - start)/(fin - start + 1)*MapW
                    coord= i,TCGA_ResultsOffset_C + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH,\
                             i + CpG_SmallMapsOverdraw,\
                             TCGA_ResultsOffset_C + 0.25*SupplCanvasH + TCGA_DrawN*SupplCanvasH + SupplCanvasH*0.3
                    IDraw.rectangle(xy = coord,fill=RGBColor)

                TextOut = '%s:'%(CrossHub_data.SamplingInfo.Dataset_Abbreviation)
                if len(CurrentCpGs) > 0:
                    TextOut += ' %d T, %d N samples (%d paired)'%(
                        CrossHub_data.SamplingInfo.Methyl_samples_count_T,
                        CrossHub_data.SamplingInfo.Methyl_samples_count_C,
                        CrossHub_data.SamplingInfo.Methyl_samples_count_paired
                    )
                if Gene in CrossHub_data.RNASeq_BasicInfo_ByGeneName:
                    if not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled) and not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired):
                        TextOut += ', RNA-Seq LogFC = %.1f / %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled,CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired)
                    elif not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled):
                        TextOut += ', RNA-Seq LogFC pooled = %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_pooled)
                    elif not math.isnan(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired):
                        TextOut += ', RNA-Seq LogFC paired = %.1f'%(CrossHub_data.RNASeq_BasicInfo_ByGeneName[Gene].LogFC_paired)
                TextOut += ' [%s]'%(CrossHub_data.SamplingInfo.Dataset_Name)
                IDraw.text(xy=(20,TCGA_ResultsOffset_C + 0.02*SupplCanvasH + TCGA_DrawN*SupplCanvasH),text=TextOut,fill=(0,0,0),font=IFont)
                # TCGA_DrawN += 1

                # рисование графика с корреляциями
                if CorrelationsPresent:
                    Cor_plot_rel_size = 0.35
                    CurrentH_Offset = TCGA_ResultsOffset_C + 0.6*SupplCanvasH + TCGA_DrawN*SupplCanvasH
                    coord = 0,CurrentH_Offset,MapW,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size
                    IDraw.rectangle(xy = coord,fill=(242,250,255))
                    coord = 0,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,MapW,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2
                    IDraw.rectangle(xy = coord,fill=(236,236,198))
                    for CpG_n in range(len(CurrentCpGs)):
                        CpG = CurrentCpGs[CpG_n]
                        add_text = None
                        if Gene in CpG.Rs_paired:
                            Rs_paired = CpG.Rs_paired[Gene]
                        else:
                            Rs_paired = CpG.Rs_paired[sorted(CpG.Rs_paired.keys())[0]]
                            add_text = 'correlations are shown for %s gene'%(sorted(CpG.Rs_paired.keys())[0])


                        # Rs_paired = (fabs(Rs_paired)**0.75)*np.sign(Rs_paired)
                        if Gene in CpG.Rs_pooled:
                            Rs_pooled = CpG.Rs_pooled[Gene]
                        else:
                            Rs_pooled = CpG.Rs_pooled[sorted(CpG.Rs_pooled.keys())[0]]
                            add_text = 'correlations are shown for %s gene'%(sorted(CpG.Rs_pooled.keys())[0])
                        # Rs_pooled = (fabs(Rs_pooled)**0.75)*np.sign(Rs_pooled)

                        i = (CurrentCpG_Positions[CpG_n] - start)/(fin - start + 1)*MapW
                        if not math.isnan(Rs_pooled):
                            RGBColor=(68,120,251)
                            coord= i ,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,i,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2 - SupplCanvasH*Cor_plot_rel_size/2*Rs_pooled
                            IDraw.rectangle(xy = coord,fill=RGBColor)

                        if not math.isnan(Rs_paired):
                            RGBColor=(241,126,56)
                            coord= i+1 ,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2,i+1,CurrentH_Offset + SupplCanvasH*Cor_plot_rel_size/2 - SupplCanvasH*Cor_plot_rel_size/2*Rs_paired
                            IDraw.rectangle(xy = coord,fill=RGBColor)

                TCGA_DrawN += 1




    # рисование профиля CG процентов
        if not opts.OnlyVisualizationOfTCGAData:
            PrevCanvasH=MainCanvasH
            CG_profileSource=list([ float(x=='C' or x=='G') for x in SourceSeq ])
            SmoothRange=7

            CG_profileSourceSmooth=smoothListGaussian(CG_profileSource,degree=SmoothRange)

            CG_profileReady=scipy.ndimage.zoom(CG_profileSourceSmooth, MapW/len(CG_profileSourceSmooth), order=4)
            CG_profileReady=list([max(min(x,1.0),0.0) for x in CG_profileReady])
            for P in range(MapW):
                if CG_profileReady[P]<0.35:
                    Div=0.35-CG_profileReady[P]
                elif CG_profileReady[P]>0.7:
                    Div=CG_profileReady[P]-0.7
                else:
                    Div=0
                MaxDiv=0.35
                RColorBase1=157
                GColorBase1=238
                BColorBase1=77
                RColorBase2=235
                GColorBase2=97
                BColorBase2=70
                RColor=int(RColorBase1+(RColorBase2-RColorBase1)*Div/MaxDiv)
                GColor=int(GColorBase1+(GColorBase2-GColorBase1)*Div/MaxDiv)
                BColor=int(BColorBase1+(BColorBase2-BColorBase1)*Div/MaxDiv)
                CColor=color(RColor,GColor,BColor)

                # Canvas.create_line(coord,fill=CColor)
                coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CG_profileReady[P]
                IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)
                # IDraw.point((P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CG_profileReady[P]),fill=(128,128,128))


            coord=0,PrevCanvasH+CG_profileH*0.9,MapW,PrevCanvasH+CG_profileH*0.9
            IDraw.line(coord,fill=(128,128,128),width=1)
            text='CG content profile (after bis., CpG meth.)'
            IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

            PrevCanvasH+=CG_profileH

        # рисование профиле CG
            if Compl:
                CG_profileSource=list([ float(x=='C' or x=='G') for x in SourceSeq.replace('CG','CA') ])
            else:
                CG_profileSource=list([ float(x=='C' or x=='G') for x in SourceSeq.replace('CG','TG') ])
            SmoothRange=7

            CG_profileSourceSmooth=smoothListGaussian(CG_profileSource,degree=SmoothRange)

            CG_profileReady=scipy.ndimage.zoom(CG_profileSourceSmooth, MapW/len(CG_profileSourceSmooth), order=4)
            CG_profileReady=list([max(min(x,1.0),0.0) for x in CG_profileReady])
            for P in range(MapW):
                if CG_profileReady[P]<0.35:
                    Div=0.35-CG_profileReady[P]
                elif CG_profileReady[P]>0.7:
                    Div=CG_profileReady[P]-0.7
                else:
                    Div=0
                MaxDiv=0.35
                RColorBase1=157
                GColorBase1=238
                BColorBase1=77
                RColorBase2=235
                GColorBase2=97
                BColorBase2=70
                RColor=int(RColorBase1+(RColorBase2-RColorBase1)*Div/MaxDiv)
                GColor=int(GColorBase1+(GColorBase2-GColorBase1)*Div/MaxDiv)
                BColor=int(BColorBase1+(BColorBase2-BColorBase1)*Div/MaxDiv)
                CColor=color(RColor,GColor,BColor)

                # Canvas.create_line(coord,fill=CColor)
                coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CG_profileReady[P]
                IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)
                # IDraw.point((P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CG_profileReady[P]),fill=(128,128,128))


            coord=0,PrevCanvasH+CG_profileH*0.9,MapW,PrevCanvasH+CG_profileH*0.9
            IDraw.line(coord,fill=(128,128,128),width=1)
            text='CG content profile (after bis., CpG ummeth.)'
            IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

            PrevCanvasH+=CG_profileH

        # рисование профиля повторностей букв
            RepeatedLettersProfile=[0.0]*len(SourceSeq)
            # SSourceSeq=list(SourceSeq)
            RepLetterCount=0
            for Pos in range(1,len(SourceSeq)):
                if SourceSeq[Pos]==SourceSeq[Pos-1]:
                    RepLetterCount+=1
                elif RepLetterCount>0:
                    for RPos in range(Pos+1-RepLetterCount,Pos+1):
                        RepeatedLettersProfile[RPos]=float(RepLetterCount)
                    RepLetterCount=1

            RepeatedLettersProfileResize=scipy.ndimage.zoom(RepeatedLettersProfile, MapW/len(RepeatedLettersProfile), order=4)
            RepeatedLettersProfileReady=list([min(x,12) for x in RepeatedLettersProfileResize])
            coord=0,PrevCanvasH+CoverageProfileH*0.9-CoverageProfileH*0.8*4/12,MapW,PrevCanvasH+CoverageProfileH*0.9-CoverageProfileH*0.8*4/12
            IDraw.line(coord,fill=(225,225,225),width=1)

            for P in range(MapW):
                CC=RepeatedLettersProfileReady[P]
                CC=min(CC,12)
                if CC<=4:
                    RColor=157
                    GColor=238
                    BColor=113
                else:
                    if CC >8:
                        DCC=8
                    else:
                        DCC=CC
                    RColorBase1=255
                    GColorBase1=231
                    BColorBase1=65
                    RColorBase2=240
                    GColorBase2=98
                    BColorBase2=83
                    RColor=max(min(int(RColorBase1+(RColorBase2-RColorBase1)*(DCC-4)/4),255),0)
                    GColor=max(min(int(GColorBase1+(GColorBase2-GColorBase1)*(DCC-4)/4),255),0)
                    BColor=max(min(int(BColorBase1+(BColorBase2-BColorBase1)*(DCC-4)/4),255),0)
                CColor=color(RColor,GColor,BColor)

                # Canvas.create_line(coord,fill=CColor)
                coord=P,PrevCanvasH+CoverageProfileH*0.9,P,PrevCanvasH+CoverageProfileH*0.9-CoverageProfileH*0.8*CC/12
                IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)
                # IDraw.point((P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CG_profileReady[P]),fill=(128,128,128))



            coord=0,PrevCanvasH+CoverageProfileH*0.9,MapW,PrevCanvasH+CoverageProfileH*0.9
            IDraw.line(coord,fill=(128,128,128),width=1)
            text='repeated letters profile; 4 let. threshold line is marked as light-grey'
            IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

            PrevCanvasH+=CG_profileH



        #  рисование профиля покрытия - ненормированные праймеры

            for DesiredXPrimers,BulkXPrimers,Type in zip([DesiredForwardPrimers,DesiredReversePrimers],[BulkForwardPrimers,BulkReversePrimers],['fwd.','rev.']):

                DesiredCoverageProfileSource=[0.0]*len(SourceSeq)
                for Pr in DesiredXPrimers:
                    for Pos in range(Pr.Start,Pr.End+1):
                        DesiredCoverageProfileSource[Pos]+=1.0

                BulkCoverageProfileSource=[0.0]*len(SourceSeq)
                for Pr in BulkXPrimers:
                    for Pos in range(Pr.Start,Pr.End+1):
                        BulkCoverageProfileSource[Pos]+=1.0

                # print(DesiredCoverageProfileSource)
                # print(BulkCoverageProfileSource)
                # exit()

                # OutFile=open('%s. bulk coverage.txt'%Type,'w')
                # Pos=0
                # for C in BulkCoverageProfileSource:
                #     OutFile.write(str(Pos)+'\t'+str(C)+'\n')
                #     Pos+=1
                # OutFile.close()

                BulkCoverageProfileResize=scipy.ndimage.zoom(BulkCoverageProfileSource, MapW/len(BulkCoverageProfileSource), order=4)
                DesiredCoverageProfileResize=scipy.ndimage.zoom(DesiredCoverageProfileSource, MapW/len(DesiredCoverageProfileSource), order=4)
                DesiredMaxV=max(max(DesiredCoverageProfileResize),1)
                BulkMaxV=max(max(BulkCoverageProfileResize),1)
                BulkCoverageProfileNorm=list([x/BulkMaxV for x in BulkCoverageProfileResize])
                DesiredCoverageProfileNorm=list([x/DesiredMaxV for x in DesiredCoverageProfileResize])
                ResizeFactor=BulkMaxV/DesiredMaxV

                # OutFile=open('%s. bulk coverage resize.txt'%Type,'w')
                # Pos=0
                # for C in BulkCoverageProfileResize:
                #     OutFile.write(str(Pos)+'\t'+str(C)+'\n')
                #     Pos+=1
                # OutFile.close()



                PrevCoord=(0,PrevCanvasH+CG_profileH*0.9)
                # print(len(BulkCoverageProfileResize))
                # print(len(DesiredCoverageProfileResize))
                for P in range(MapW):
                    BulkCC=min(BulkCoverageProfileNorm[P],0.7)
                    RColorBase1=255
                    GColorBase1=231
                    BColorBase1=65
                    RColorBase2=104
                    GColorBase2=204
                    BColorBase2=250
                    RColor=int(RColorBase1+(RColorBase2-RColorBase1)*BulkCC/0.7)
                    GColor=int(GColorBase1+(GColorBase2-GColorBase1)*BulkCC/0.7)
                    BColor=int(BColorBase1+(BColorBase2-BColorBase1)*BulkCC/0.7)

                    BulkRColor=int((255+RColor*1.5)/2.5)
                    BulkGColor=int((255+GColor*1.5)/2.5)
                    BulkBColor=int((255+BColor*1.5)/2.5)

                    DesiredCC=min(DesiredCoverageProfileNorm[P],0.7)
                    RColorBase1=255
                    GColorBase1=231
                    BColorBase1=65
                    RColorBase2=104
                    GColorBase2=204
                    BColorBase2=250
                    DesiredRColor=int(RColorBase1+(RColorBase2-RColorBase1)*BulkCC/0.7)
                    DesiredGColor=int(GColorBase1+(GColorBase2-GColorBase1)*BulkCC/0.7)
                    DesiredBColor=int(BColorBase1+(BColorBase2-BColorBase1)*BulkCC/0.7)


                    coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*BulkCoverageProfileNorm[P]
                    IDraw.line(coord,fill=(BulkRColor,BulkGColor,BulkBColor),width=1)

                    coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*DesiredCoverageProfileNorm[P]/ResizeFactor
                    IDraw.line(coord,fill=(DesiredRColor,DesiredGColor,DesiredBColor),width=1)

                    NewCoord=(P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*DesiredCoverageProfileNorm[P]/ResizeFactor)
                    IDraw.line(PrevCoord+NewCoord,fill=(128,128,128),width=1)
                    PrevCoord=NewCoord


                coord=0,PrevCanvasH+CoverageProfileH*0.9,MapW,PrevCanvasH+CoverageProfileH*0.9
                IDraw.line(coord,fill=(128,128,128),width=1)
                text='%s primers coverage profile: with desired param. (dark fill); with limit param. (light fill), max. %d'%(Type,int(BulkMaxV))
                IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

                PrevCanvasH+=CG_profileH


        #  рисование профиля покрытия - нормализованные праймеры
            for XPrimers,Type in zip([ForwardPrimers,ReversePrimers],['fwd.','rev.']):
                CoverageProfileSource=[0.0]*len(SourceSeq)
                for Pr in XPrimers:
                    for Pos in range(Pr.Start,Pr.End+1):
                        CoverageProfileSource[Pos]+=1.0

                CoverageProfileResize=scipy.ndimage.zoom(CoverageProfileSource, MapW/len(CoverageProfileSource), order=4)
                MaxV=max(max(CoverageProfileResize),1)
                CoverageProfileNorm=list([x/MaxV for x in CoverageProfileResize])

                # OutFile=open('%s. normal. coverage.txt'%Type,'w')
                # Pos=0
                # for C in CoverageProfileSource:
                #     OutFile.write(str(Pos)+'\t'+str(C)+'\n')
                #     Pos+=1
                # OutFile.close()
                #
                # OutFile=open('%s. normal. coverage resize.txt'%Type,'w')
                # Pos=0
                # for C in CoverageProfileResize:
                #     OutFile.write(str(Pos)+'\t'+str(C)+'\n')
                #     Pos+=1
                # OutFile.close()


                # print(len(CoverageProfileResize))
                # PrevCoord=(0,PrevCanvasH+CG_profileH*0.9)
                for P in range(MapW):
                    CC=min(CoverageProfileNorm[P],0.8)
                    RColorBase1=255
                    GColorBase1=231
                    BColorBase1=65
                    RColorBase2=104
                    GColorBase2=204
                    BColorBase2=250
                    RColor=int(RColorBase1+(RColorBase2-RColorBase1)*CC/0.8)
                    GColor=int(GColorBase1+(GColorBase2-GColorBase1)*CC/0.8)
                    BColor=int(BColorBase1+(BColorBase2-BColorBase1)*CC/0.8)

                    coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*CoverageProfileNorm[P]
                    IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)

                    # if CoverageProfileNorm[P] < 0.05:
                    #     print('Type=%s; Cov prof norm=%.3f, Pos=%d'%(Type,CoverageProfileNorm[P],P))
                    #     exit()


                coord=0,PrevCanvasH+CoverageProfileH*0.9,MapW,PrevCanvasH+CoverageProfileH*0.9
                IDraw.line(coord,fill=(128,128,128),width=1)
                text='normalized %s primers coverage profile, max. %d'%(Type,int(MaxV))
                IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

                PrevCanvasH+=CG_profileH

        #  рисование профиля штрафов с фоновым отображением областей, вошедших в штрафы

            if not opts.BypassSpecificityTest:
                CoverageProfileSource=[0.0]*len(SourceSeq)
                for Pr in SpecTestForwardPrimers + SpecTestReversePrimers:
                    for Pos in range(Pr.Start,Pr.End+1):
                        CoverageProfileSource[Pos]+=1.0


                PenaltyProfileSource=[0.0]*len(SourceSeq)
                for Pr in SpecTestForwardPrimers + SpecTestReversePrimers:
                    for Pos in range(Pr.Start,Pr.End+1):
                        PenaltyProfileSource[Pos] += Pr.SpecPenalty

                for Pos in range(len(SourceSeq)):
                    if CoverageProfileSource[Pos] > 0:
                        PenaltyProfileSource[Pos] = (PenaltyProfileSource[Pos]/CoverageProfileSource[Pos])**0.5
                    else:
                        PenaltyProfileSource[Pos] = 0

                CoverageProfileResize=scipy.ndimage.zoom(CoverageProfileSource, MapW/len(CoverageProfileSource), order=4)

                PenaltyProfileResize=scipy.ndimage.zoom(PenaltyProfileSource, MapW/len(PenaltyProfileSource), order=4)
                MaxV=max(max(PenaltyProfileResize),1)
                PenaltyProfileNorm=list([x/MaxV for x in PenaltyProfileResize])

                for P in range(MapW):
                    coord=P,PrevCanvasH+CoverageProfileH*0.9,P,PrevCanvasH+CoverageProfileH*0.1
                    if CoverageProfileResize[P] > 0:
                        IDraw.line(coord,fill=(242,253,235),width=1)
                    else:
                        IDraw.line(coord,fill=(240,210,192),width=1)

                PrevCoord=None
                for P in range(MapW):
                    CC=min(PenaltyProfileNorm[P],0.8)
                    RColorBase1=196
                    GColorBase1=235
                    BColorBase1=148
                    RColorBase2=255
                    GColorBase2=182
                    BColorBase2=133
                    RColor=int(RColorBase1+(RColorBase2-RColorBase1)*CC/0.8)
                    GColor=int(GColorBase1+(GColorBase2-GColorBase1)*CC/0.8)
                    BColor=int(BColorBase1+(BColorBase2-BColorBase1)*CC/0.8)

                    if CoverageProfileResize[P]>0:
                        coord=P,PrevCanvasH+CG_profileH*0.9,P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*PenaltyProfileNorm[P]
                        IDraw.line(coord,fill=(RColor,GColor,BColor),width=1)
                        if PrevCoord==None:
                            PrevCoord=P,PrevCanvasH+CG_profileH*0.9-CG_profileH*0.8*PenaltyProfileNorm[P]
                            IDraw.line(coord,fill=(128,128,128),width=1)
                        else:
                            IDraw.line(PrevCoord+coord[2:],fill=(128,128,128),width=1)
                            PrevCoord = coord[2:]
                    else:
                        if PrevCoord!=None:
                            IDraw.line(PrevCoord+(P-1,PrevCanvasH+CG_profileH*0.9),fill=(128,128,128),width=1)
                        PrevCoord=None




                text='specificity test penalty profile (sqrt, max %d). non-tested areas have reddish background'%(int(MaxV**2))
                IDraw.text(xy=(20, PrevCanvasH),text=text,fill=(0,0,0),font=IFont)

                PrevCanvasH+=CG_profileH

        # рисование непосредственно ампликонов

            MaxP=200 # Максимального штрафа    Цвет
            for PN in range(len(CompletePairsList)):
                P=CompletePairsList[PN]
                SupplCanvasHStart=PrevCanvasH+PN*SupplCanvasH
                SupplCanvasWStart=MainCanvasW/(InChrEnd-InChrStart)*(P.ForwardPrimer.InChrStart-InChrStart)
                SupplCanvasW=MainCanvasW/(InChrEnd-InChrStart)*(P.ReversePrimer.InChrEnd - P.ForwardPrimer.InChrStart)

            # вычисление основного цвета для ампликона
                PCount = min(P.Penalty,MaxP)
                RColorBase1=183
                GColorBase1=234
                BColorBase1=124
                RColorBase2=237
                GColorBase2=102
                BColorBase2=49
                RColorBack=int(RColorBase1+(RColorBase2-RColorBase1)*PCount/MaxP)
                GColorBack=int(GColorBase1+(GColorBase2-GColorBase1)*PCount/MaxP)
                BColorBack=int(BColorBase1+(BColorBase2-BColorBase1)*PCount/MaxP)
                PColor=color(RColorBack,GColorBack,BColorBack)

                coord=SupplCanvasWStart,SupplCanvasHStart+SupplCanvasH*0.1,SupplCanvasWStart+SupplCanvasW,SupplCanvasHStart+SupplCanvasH*0.45
                # Canvas.create_rectangle(coord,fill=PColor, width=0)
                IDraw.rectangle(xy=coord,fill=(RColorBack,GColorBack,BColorBack))

                # добавление полосок с ЦпГ

                MaxS = 1.5
                for i in range(int(SupplCanvasWStart),int(SupplCanvasWStart+SupplCanvasW)):
                    # if SmallMap[i]>0:
                    SCount = min(SmallMap[i],MaxS)
                    RColorBase1=RColorBack
                    GColorBase1=GColorBack
                    BColorBase1=BColorBack
                    RColorBase2=0
                    GColorBase2=126
                    BColorBase2=255
                    RColor=int(RColorBase1+float(RColorBase2-RColorBase1)*float(SCount)/float(MaxS))
                    GColor=int(GColorBase1+float(GColorBase2-GColorBase1)*float(SCount)/float(MaxS))
                    BColor=int(BColorBase1+float(BColorBase2-BColorBase1)*float(SCount)/float(MaxS))
                    CColor=color(RColor,GColor,BColor)
                    coord=i,SupplCanvasHStart+SupplCanvasH*0.1,i,SupplCanvasHStart+SupplCanvasH*0.45

                    # Canvas.create_line(coord,fill=CColor)
                    IDraw.line(xy=coord,fill=(RColor,GColor,BColor),width=1)

                #  рисую пряаймеры

                PrimerWStart=MainCanvasW/(InChrEnd-InChrStart)*(P.ForwardPrimer.InChrStart-InChrStart)
                PrimerWEnd=MainCanvasW/(InChrEnd-InChrStart)*(P.ForwardPrimer.InChrEnd-InChrStart)
                coord=PrimerWStart,SupplCanvasHStart+SupplCanvasH*0.02,PrimerWEnd,SupplCanvasHStart+SupplCanvasH*0.08
                # Canvas.create_rectangle(coord,fill='#9b9b9b', width=0)
                IDraw.rectangle(xy=coord,fill=(155,155,155))

                PrimerWStart=MainCanvasW/(InChrEnd-InChrStart)*(P.ReversePrimer.InChrStart-InChrStart)
                PrimerWEnd=MainCanvasW/(InChrEnd-InChrStart)*(P.ReversePrimer.InChrEnd-InChrStart)
                coord=PrimerWStart,SupplCanvasHStart+SupplCanvasH*0.47,PrimerWEnd,SupplCanvasHStart+SupplCanvasH*0.53
                # Canvas.create_rectangle(coord,fill='#9b9b9b', width=0)
                IDraw.rectangle(xy=coord,fill=(155,155,155))

                AmpLength=P.ReversePrimer.InChrEnd-P.ForwardPrimer.InChrStart+1

                text='id.%d, amp.%d, p.%.1f, CpG in pr.: %d/%d'%(PN+1,AmpLength,P.Penalty,P.ForwardPrimer.CpGsCount,P.ReversePrimer.CpGsCount)
                w,h = IDraw.textsize(text=text,font=IFont)
                # print(w)
                # print(SupplCanvasWStart+5)
                # print(MapW)
                IDraw.text(xy=(min(SupplCanvasWStart+5,MapW-5-w), SupplCanvasHStart+SupplCanvasH*0.57),text=text,fill=(0,0,0),font=IFont)

                # print(P.ReversePrimer.Sequence)
                # L=len(SourceSeq)
                # print(SourceSeq[P.ReversePrimer.Start:P.ReversePrimer.End+1])
                # print(SourceSeq[L-P.ReversePrimer.End:L-P.ReversePrimer.Start+1])
                # Canvas.create_text(SupplCanvasWStart+5, SupplCanvasHStart+SupplCanvasH*0.57, anchor=NW,font="Purisa",text=text)

        if opts.CustomResultsName == None:
            if GeneCoordinates.IsFragmentOfLargeGene:
                if not os.path.exists(os.path.join('Run_%s.files'%RunTime,'Fragments')):
                    os.mkdir(os.path.join('Run_%s.files'%RunTime,'Fragments'))
                ImageFileName = os.path.join('Run_%s.files'%RunTime,'Fragments','id.%d_%s__Chr.%s___%d-%d.png'%(CurrentGeneID+1,Gene,Chr,start,fin))
            else:
                ImageFileName = os.path.join('Run_%s.files'%RunTime,'id.%d_%s__Chr.%s___%d-%d.png'%(CurrentGeneID+1,Gene,Chr,start,fin))
        else:
            if GeneCoordinates.IsFragmentOfLargeGene:
                if not os.path.exists(os.path.join('%s.files'%opts.CustomResultsName,'Fragments')):
                    os.mkdir(os.path.join('%s.files'%opts.CustomResultsName,'Fragments'))
                ImageFileName = os.path.join('%s.files'%opts.CustomResultsName,'Fragments','id.%d_%s__Chr.%s___%d-%d.png'%(CurrentGeneID+1,Gene,Chr,start,fin))
            else:
                ImageFileName = os.path.join('%s.files'%opts.CustomResultsName,'id.%d_%s__Chr.%s___%d-%d.png'%(CurrentGeneID+1,Gene,Chr,start,fin))
        # Canvas.pack()
        # Canvas.postscript(file='test1.ps', colormode='color')
        del IDraw
        IImage.save(ImageFileName,format='PNG')
        # image1 = Image.new("RGB", (MainCanvasW, MainCanvasH+SupplCanvasH*len(CompletePairsList)), (255,255,255))
        # image1.save('TEST1.jpg')

        # Canvas.fil

        # top.mainloop()

        PrimerInfoFile.write("<font size=6.4><i>%s</i></font><br>"%Gene)
        if ImageFileName!=None:
            PrimerInfoFile.write("<img src=\"%s\"><br><br>" %ImageFileName)
        if Compl:
            TextInsert=' (complement strand)'
        else:
            TextInsert=''

        if(len(SegDescriptions)>0):
            PrimerInfoFile.write('<b>ENCODE Genome segments color legend:</b><br>\n')
            PrimerInfoFile.write('<table>')
            PrimerInfoFile.write('<tr><td>Active Promoter</td><td><span style="background-color: #ff0000">Bright red</span></td><td>Predicted promoter region including TSS</td></tr>\n') #ff0000'
            PrimerInfoFile.write('<tr><td>Promoter Flanking</td><td><span style="background-color: #f98080">Light red</span></td><td>Predicted promoter flanking region</td></tr>\n') #f98080'
            PrimerInfoFile.write('<tr><td>Inactive/Poised Promoter</td><td><span style="background-color: #e89dff">Purple</span></td><td>Predicted inactive/poised promoter</td></tr>\n') #e89dff'
            PrimerInfoFile.write('<tr><td>Strong Enhancer</td><td><span style="background-color: #faca00">Orange</span></td><td>Predicted enhancer</td></tr>\n') #faca00'
            PrimerInfoFile.write('<tr><td>Weak Enhancer</td><td><span style="background-color: #fff29b">Light yellow</span></td><td>Predicted weak enhancer or open chromatin cis regulatory element</td></tr>\n') #ffe641'
            PrimerInfoFile.write('<tr><td>Insulator</td><td><span style="background-color: #82aafa">Blue</span></td><td>Predicted insulator (CTCF enriched region)</td></tr>\n') #82aafa'
            PrimerInfoFile.write('<tr><td>Transcription associated region</td><td><span style="background-color: #01b954">Deep green</span></td><td>Predicted transcribed region</td></tr>\n') #01b954'
            PrimerInfoFile.write('<tr><td>Low Activity</td><td><span style="background-color: #a8fca8">Light green</span></td><td>Predicted low activity region</td></tr>\n') #a8fca8'
            PrimerInfoFile.write('<tr><td>Polycomb Repressed</td><td><span style="background-color: #9d9d9d">Gray</span></td><td>Predicted repressed region</td></tr>\n') #9d9d9d'
            PrimerInfoFile.write('<tr><td>Heterochromatin/unknown</td><td><span style="background-color: #f5f5f5">White</span></td><td>Heterochromatin/unknown/ </td></tr>\n') #f5f5f5'
            PrimerInfoFile.write('</table><br><br>\n')


        # PrimerInfoFile.write("<font size=4.5><a href=\"http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%s%%3A%d-%d&hgsid=199006321_NVApKptlIQZyk3erAPDwgUyr7MGv\">go to UCSC genome browser%s</a></font><br><br>"%(Chr,InChrStart,InChrEnd,TextInsert))
        PrimerInfoFile.write("<font size=4.5><a href=\"http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=%s&position=chr%s%%3A%d-%d\">go to UCSC genome browser%s</a></font><br><br>"%(opts.GenomeVersion,Chr,InChrStart,InChrEnd,TextInsert))
        if not opts.OnlyVisualizationOfTCGAData:
            for PN in range(len(CompletePairsList)):
                PrimerInfoFile.write("<br><hr style=\"border: 1px dotted #2c1294;\"><br>")
                PrimerInfoFile.write("<font size=5.3><b><i>%s</i>, fragment %d"%(Gene, PN+1))
                PrimerInfoFile.write("</font></b><br><br>")
                PrimerInfoFile.write("<table border='1'>")
                PrimerInfoFile.write("<b><tr><td></td><td>Sequence</td><td>Length</td><td>Start</td><td>End</td><td>Tm</td><td>comments</td><td>dimers</td></tr></b>")
                P=CompletePairsList[PN]
                Init=P.ForwardPrimer.Sequence
                Mod=P.ForwardPrimer.Sequence.replace('CG','TG')

                FPrimerSeq=P.ForwardPrimer.Sequence
                if opts.Mode == 'Methyl':  FPrimerSeqUnmeth = P.ForwardPrimer.SequenceUnmeth
                if opts.Mode == 'Simple':  FPrimerSeqUnmeth = P.ForwardPrimer.Sequence
                RPrimerSeq=P.ReversePrimer.Sequence
                if opts.Mode == 'Methyl':  RPrimerSeqUnmeth = P.ReversePrimer.SequenceUnmeth
                if opts.Mode == 'Simple':  RPrimerSeqUnmeth = P.ReversePrimer.Sequence
                AmpSize=P.ReversePrimer.InChrEnd - P.ForwardPrimer.InChrStart + 1

                DimersC=MaximumEnergy([FPrimerSeq,FPrimerSeqUnmeth],GenPar.TDimersdGPenaltyIncrement,GenPar.FDimersdGPenaltyIncrement)
                PrimerCpG_Positions = []
                Pos=-1
                while True:
                    Pos = P.ForwardPrimer.SequencePlus.find('CG',Pos+1)
                    if Pos >= 0:
                        PrimerCpG_Positions.append(Pos)
                        PrimerCpG_Positions.append(Pos+1)
                    else:
                        break

                PrimerCpG_Positions = [x - P.ForwardPrimer.SeqStartInSeqPlus for x in PrimerCpG_Positions]

                HTML_PrimerSeqs=[]
                for CSeq in list(sorted(set([FPrimerSeq,FPrimerSeqUnmeth]))):
                    LettersList = list(CSeq)
                    for N in range(len(LettersList)):
                        if N in PrimerCpG_Positions:
                            LettersList[N] = '<span style="background-color: #BBCCFF">%s</span>'%LettersList[N]
                    HTML_PrimerSeqs.append(''.join(LettersList))


                PrimerInfoFile.write("<tr><td>forward</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%d CpG</td><td>%s<br>%s<br>%s<br>dG = %.1f kc/m</td></tr>\n"
                    %('<br>'.join(HTML_PrimerSeqs),
                        P.ForwardPrimer.InChrEnd-P.ForwardPrimer.InChrStart+1,
                        intWithCommas(P.ForwardPrimer.InChrStart),
                        intWithCommas(P.ForwardPrimer.InChrEnd),
                        '<br>'.join(['%.1f C'%ThermalFunction(x) for x in sorted((set([FPrimerSeq,FPrimerSeqUnmeth]   )))  ]),
                        P.ForwardPrimer.CpGsCount,
                        DimersC.upper.replace(' ','&nbsp'),
                        DimersC.binds.replace(' ','&nbsp'),
                        DimersC.lower.replace(' ','&nbsp'),
                        DimersC.Energy))

                DimersC=MaximumEnergy([RPrimerSeq,RPrimerSeqUnmeth],GenPar.TDimersdGPenaltyIncrement,GenPar.FDimersdGPenaltyIncrement)
                PrimerCpG_Positions = []
                Pos=-1
                while True:
                    Pos = P.ReversePrimer.SequencePlus.find('CG',Pos+1)
                    if Pos >= 0:
                        PrimerCpG_Positions.append(Pos)
                        PrimerCpG_Positions.append(Pos+1)
                    else:
                        break

                PrimerCpG_Positions = [x - P.ReversePrimer.SeqStartInSeqPlus for x in PrimerCpG_Positions]

                HTML_PrimerSeqs=[]
                for CSeq in list(sorted(set([RPrimerSeq,RPrimerSeqUnmeth]))):
                    LettersList = list(CSeq)
                    for N in range(len(LettersList)):
                        if N in PrimerCpG_Positions:
                            LettersList[N] = '<span style="background-color: #BBCCFF">%s</span>'%LettersList[N]
                    HTML_PrimerSeqs.append(''.join(LettersList))



                PrimerInfoFile.write("<tr><td>reverse</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%d CpG</td><td>%s<br>%s<br>%s<br>dG = %.1f kc/m</td></tr>\n"
                    %('<br>'.join(HTML_PrimerSeqs),
                        P.ReversePrimer.InChrEnd-P.ReversePrimer.InChrStart+1,
                        intWithCommas(P.ReversePrimer.InChrStart),
                        intWithCommas(P.ReversePrimer.InChrEnd),
                        '<br>'.join(['%.1f C'%ThermalFunction(x) for x in sorted((set(  [RPrimerSeq,RPrimerSeqUnmeth]   )))  ]),
                        P.ReversePrimer.CpGsCount,
                        DimersC.upper.replace(' ','&nbsp'),
                        DimersC.binds.replace(' ','&nbsp'),
                        DimersC.lower.replace(' ','&nbsp'),
                        DimersC.Energy))
                PrimerInfoFile.write("</table><br><table border='1'>\n")
                PrimerInfoFile.write("<b><tr><td>Amplicon size</td><td>Amplicon Tm</td><td>CpG count</td><td>Penalty</td><td>False PCR products:<br>sizes and mismatches in primers</td><td>Primer max.en. dimer</td></tr> </b>")

                DimersC=MaximumEnergy([FPrimerSeq,FPrimerSeqUnmeth,RPrimerSeq,RPrimerSeqUnmeth],GenPar.TDimersdGPenaltyIncrement,GenPar.FDimersdGPenaltyIncrement)
                PrimerInfoFile.write("<tr><td>%d bp.</td><td>&nbsp%.1f C&nbsp</td><td>%d CpG</td><td>%.1f</td>"%(AmpSize,ThermalFunction(P.AmpliconSequence),P.CpGsCount,P.Penalty))
                if P.ObjectiveSpecTest:
                    if P.CompleteFalseAmpliconsCount==0:
                        PrimerInfoFile.write('<td><i>None</i></td>')
                    elif P.CompleteFalseAmpliconsCount <= GenPar.FalseAmpliconsToShow:
                        PrimerInfoFile.write('<td>total %d;<br>%s</td>'%(P.CompleteFalseAmpliconsCount,'<br>'.join(['%d bp: %d/%d'%(x.Size,x.ForwardPrimerMismatches,x.ReversePrimerMismatches ) for x in P.FalseAmplicons])))
                    else:
                        PrimerInfoFile.write('<td><i>total %d; displaying first %d</i><br>%s</td>'%(P.CompleteFalseAmpliconsCount,GenPar.FalseAmpliconsToShow,'<br>'.join(['%d bp: %d/%d'%(x.Size,x.ForwardPrimerMismatches,x.ReversePrimerMismatches ) for x in P.FalseAmplicons[:GenPar.FalseAmpliconsToShow]])))
                else:
                    PrimerInfoFile.write('<td>%s</td>'%P.NonObjectiveSpecTestInfo)
                PrimerInfoFile.write('<td>%s<br>%s<br>%s<br>dG = %1.f kc/m</td></tr>\n'%(DimersC.upper.replace(' ','&nbsp'),
                                        DimersC.binds.replace(' ','&nbsp'),
                                        DimersC.lower.replace(' ','&nbsp'),
                                        DimersC.Energy))
                PrimerInfoFile.write("</table><br><br>")

                amplicon_wrap_strings = '<br>'.join(textwrap.wrap(P.AmpliconSequence,80))
                amplicon_wrap_strings = amplicon_wrap_strings.replace('CG','<span style="background-color: #CCDDFF">CG</span>').replace('C<br>G','<span style="background-color: #BBCCFF">C<br>G</span>')
                CG_in_amplicon_count = P.AmpliconSequence.count('CG')
                PrimerInfoFile.write("Amplicon sequence (%d CpG sites):\n<br>"%CG_in_amplicon_count + amplicon_wrap_strings + '<br>')

            PrimerInfoFile.write('<hr style="border: 2px dotted #2c1294;"><br><br><br><br>')

        CurrentGeneID+=1
        del ForwardPrimers[:]
        del ReversePrimers[:]
        del BulkForwardPrimers[:]
        del BulkReversePrimers[:]
        del SpecTestForwardPrimers[:]
        del SpecTestReversePrimers[:]
        del Pairs[:]
        del PrimersByID

        gc.collect()

    if PrimerInfoFile!=None:
        PrimerInfoFile.write("</font></body>\n</html>")

    print('Completed.')


if __name__ == "__main__":
    main()

