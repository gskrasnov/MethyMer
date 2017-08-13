__author__ = 'Administrator'
from math import log,exp


def LoadThermoParameters():
    global hAA,hAT,hAG,hAC,hTA,hTT,hTG,hTC,hGA,hGT,hGG,hGC,hCA,hCT,hCG,hCC,hI,hIgc,hIat,hIsym,hIterm
    global sAA,sAT,sAG,sAC,sTA,sTT,sTG,sTC,sGA,sGT,sGG,sGC,sCA,sCT,sCG,sCC,sI,sIgc,sIat,sIsym,sIterm
    hAA=-9100
    hAT=-8600
    hAG=-7800
    hAC=-6500
    hTA=-6000
    hTT=-9100
    hTG=-5800
    hTC=-5600
    hGA=-5600
    hGT=-6500
    hGG=-11000
    hGC=-11100
    hCA=-5800
    hCT=-7800
    hCG=-11900
    hCC=-11000
    hI=-350

    sAA=-24.0
    sAT=-23.9
    sAG=-20.8
    sAC=-17.3
    sTA=-16.9
    sTT=-24.0
    sTG=-12.9
    sTC=-13.5
    sGA=-13.5
    sGT=-17.3
    sGG=-26.6
    sGC=-26.7
    sCA=-12.9
    sCT=-20.8
    sCG=-27.8
    sCC=-26.6
    sI=-19.4


def Entalpy(primer):
    global hAA,hAT,hAG,hAC,hTA,hTT,hTG,hTC,hGA,hGT,hGG,hGC,hCA,hCT,hCG,hCC,hI,hIgc,hIat,hIsym,hIterm
    global sAA,sAT,sAG,sAC,sTA,sTT,sTG,sTC,sGA,sGT,sGG,sGC,sCA,sCT,sCG,sCC,sI,sIgc,sIat,sIsym,sIterm
    global dimers
    length=len(primer)
    H=0
    for i in range(0,length-1):
        if primer[i]=='A' and primer[i+1]=='A':
                H=H+hAA
        if primer[i]=='A' and primer[i+1]=='T':
                H=H+hAT
        if primer[i]=='A' and primer[i+1]=='G':
                H=H+hAG
        if primer[i]=='A' and primer[i+1]=='C':
                H=H+hAC
        if primer[i]=='T' and primer[i+1]=='A':
                H=H+hTA
        if primer[i]=='T' and primer[i+1]=='T':
                H=H+hTT
        if primer[i]=='T' and primer[i+1]=='G':
                H=H+hTG
        if primer[i]=='T' and primer[i+1]=='C':
                H=H+hTC
        if primer[i]=='G' and primer[i+1]=='A':
                H=H+hGA
        if primer[i]=='G' and primer[i+1]=='T':
                H=H+hGT
        if primer[i]=='G' and primer[i+1]=='G':
                H=H+hGG
        if primer[i]=='G' and primer[i+1]=='C':
                H=H+hGC
        if primer[i]=='C' and primer[i+1]=='A':
                H=H+hCA
        if primer[i]=='C' and primer[i+1]=='T':
                H=H+hCT
        if primer[i]=='C' and primer[i+1]=='G':
                H=H+hCG
        if primer[i]=='C' and primer[i+1]=='C':
                H=H+hCC
    return H

def Entropy(primer):
    S=0
    for i in range (0,len(primer)-1):
        if primer[i]=='A' and primer[i+1]=='A':
                S=S+sAA

        if primer[i]=='A' and primer[i+1]=='T':
                S=S+sAT

        if primer[i]=='A' and primer[i+1]=='G':
                S=S+sAG

        if primer[i]=='A' and primer[i+1]=='C':
                S=S+sAC

        if primer[i]=='T' and primer[i+1]=='A':
                S=S+sTA

        if primer[i]=='T' and primer[i+1]=='T':
                S=S+sTT

        if primer[i]=='T' and primer[i+1]=='G':
                S=S+sTG

        if primer[i]=='T' and primer[i+1]=='C':
                S=S+sTC

        if primer[i]=='G' and primer[i+1]=='A':
                S=S+sGA

        if primer[i]=='G' and primer[i+1]=='T':
                S=S+sGT

        if primer[i]=='G' and primer[i+1]=='G':
                S=S+sGG

        if primer[i]=='G' and primer[i+1]=='C':
                S=S+sGC

        if primer[i]=='C' and primer[i+1]=='A':
                S=S+sCA

        if primer[i]=='C' and primer[i+1]=='T':
                S=S+sCT

        if primer[i]=='C' and primer[i+1]=='G':
                S=S+sCG

        if primer[i]=='C' and primer[i+1]=='C':
                S=S+sCC

        return S




class TDimerInfo():
    def __init__(self):
        self.upper=[' ']*500
        self.lower=[' ']*500
        self.binds=[' ']*500
        self.Energy=0

def MaximumEnergy(PrimerList,TPen,FPen):
    BestDimersInfo=None
    BestEnergy=None
    for A in range(len(PrimerList)):
        for B in range(len(PrimerList)):
            CDimers=Energy(PrimerList[A],PrimerList[B],TPen,FPen)
            if CDimers==None:
                continue
            if BestEnergy==None:
                BestEnergy=CDimers.Energy
                BestDimersInfo=CDimers
            elif BestEnergy > CDimers.Energy:
                BestEnergy = CDimers.Energy
                BestDimersInfo=CDimers

    if BestDimersInfo==None:
        BestDimersInfo=TDimerInfo()
        BestDimersInfo.upper=' '
        BestDimersInfo.lower=' '
        BestDimersInfo.binds=' '
    return BestDimersInfo


def Energy(PrimerFwdInit, PrimerRewInit, TPen, FPen):
    global ThermoMethodics,SaltSwitcher,DimersTemperature,WriteSwitcher,upper,lower,binds

    upper=[' ']*500
    lower=[' ']*500
    binds=[' ']*500


    PrimerFwd=[' ']*500
    PrimerRew=[' ']*500
    for P in range(0,len(PrimerFwdInit)):
        PrimerFwd[P]=PrimerFwdInit[P]
    for P in range(0,len(PrimerRewInit)):
        PrimerRew[P]=''.join(list(reversed(PrimerRewInit)))[P]

    FwdLen=len(PrimerFwdInit)
    RewLen=len(PrimerRewInit)

    LinkLen=[0]*500
    LinkStart=[0]*500
    G=[0]*500

    EnT=0
    MaxEn=-1000000
    LinkExist=False
    for a in range (0,FwdLen):
        b=0
        LinkNumber=0

        while True:
            Condition=(((PrimerFwd[a+b]=='A') and (PrimerRew[b]=='T'))  or ((PrimerFwd[a+b]=='T') and (PrimerRew[b]=='A'))  or ((PrimerFwd[a+b]=='G') and (PrimerRew[b]=='C'))  or ((PrimerFwd[a+b]=='C') and (PrimerRew[b]=='G')))
            if Condition:
                LinkOpen=True
                LinkExist=True
                LinkLen[LinkNumber]=1
                LinkStart[LinkNumber]=a+b
                b+=1
                # if a+b<FwdLen:
                while(((PrimerFwd[a+b]=='A') and (PrimerRew[b]=='T'))  or ((PrimerFwd[a+b]=='T') and (PrimerRew[b]=='A'))  or ((PrimerFwd[a+b]=='G') and (PrimerRew[b]=='C'))  or ((PrimerFwd[a+b]=='C') and (PrimerRew[b]=='G'))):
                    LinkLen[LinkNumber]+=1
                    b+=1
                    # if a+b>=FwdLen:
                    #     break

                LinkNumber+=1
                LinkOpen=False

            else:
                LinkOpen=False
                b+=1
                if b>=min(RewLen,FwdLen-a):
                    break


        if LinkNumber==0:
            continue

        for c in range(0,LinkNumber):
            if LinkLen[c]>1:
                filetext=PrimerFwd[LinkStart[c] : LinkStart[c]+LinkLen[c]]

                H=Entalpy(filetext)
                S=Entropy(filetext)
                G[c]=H-DimersTemperature*S
            else:
                if ((PrimerFwd[LinkStart[c]]=='A')  or (PrimerFwd[LinkStart[c]]=='T')):
                    H=-3717
                    S=-10.65
                    G[c]=(H-DimersTemperature*S)/2

                else:
                    H=-4733
                    S=-11.91
                    G[c]=(H-DimersTemperature*S)/2

            if SaltSwitcher:
                G[c]=G[c]+16.6*S*log(SaltConc,2.71828)/exp(1)



        RMid=0
        e=0
        LinkMidSize=0
        for c in range(0,LinkNumber):
            for d in range(0,LinkLen[c]):
                RMid=RMid+float(LinkStart[c]+d)
                e+=1

            LinkMidSize=LinkMidSize+float(LinkLen[c])


        RMid=RMid/e
        LinkMidSize=LinkMidSize/LinkNumber
        e=0

        RModMid=0
        for c in range (0,LinkNumber):
            for d in range(0,LinkLen[c]):
                RModMid=RModMid+max(RMid,float(LinkStart[c]+d))-min(RMid,float(LinkStart[c]+d))
                e+=1


        RModMid=RModMid/e
        if(RModMid<1):
            RModMid=1

        if RModMid!=0:
            Expo=exp(log(RModMid/(float(max(LinkMidSize-1,float(1)))),2.71828)*0.75)*1.5
        else:
            Expo=RModMid

        En=0
        for c in range(0,LinkNumber):
            if G[c]<0:
                En=En+exp(Expo*log((-1)*G[c],2.71828))
            else:
                En=En-exp(Expo*log(G[c],2.71828))


        if En>0:
            En=(-1)*exp(log(En,2.71828)/Expo)

        else:
            En=exp(log((-1)*En,2.71828)/Expo)


        if a+RewLen>FwdLen:
            EnT=En*(TPen/2)*(float(RewLen-(FwdLen-a)))/(max((float(RewLen)),float(5)))

        EnF=En*(TPen/2)*(float(a))/(max(float(FwdLen),float(5)))
        En=En+EnT+EnF

                # sprintf(upper,"A:%d LinkNum: %d; Energy: %g Expo: %g RModMid: %g LinkMidSize: %g, G1: %g G2: %g",a,LinkNumber,En,Expo,RModMid,LinkMidSize,G[0],G[1]);
                # MessageBox(NULL,upper,upper,MB_OK);

        if(((-1)*En)>MaxEn):
            MaxEn=(-1)*En
            if(WriteSwitcher):
                for b in range(0,FwdLen):
                    upper[b]=PrimerFwd[b]

                upperLen=FwdLen
                for b in range(0,a+1):
                    lower[b]=' '

                for b in range(0,RewLen):
                    lower[a+b]=PrimerRew[b]

                lowerLen=a+RewLen

                for b in range (0,FwdLen):
                    if(((upper[b]=='A') and (lower[b]=='T'))  or ((upper[b]=='T') and (lower[b]=='A'))  or ((upper[b]=='G') and (lower[b]=='C'))  or ((upper[b]=='C') and (lower[b]=='G'))):
                            binds[b]='|'
                    else:
                        binds[b]=' '

                bindsLen=min(FwdLen,a+RewLen)


    PrimerFwd=[' ']*500
    for a in range(0,RewLen):
        PrimerFwd[a]=PrimerRewInit[RewLen-a-1]


    PrimerRew=[' ']*500
    for a in range(0,FwdLen):
        PrimerRew[a]=PrimerFwdInit[a]


    a=RewLen
    RewLen=FwdLen
    FwdLen=a

    LinkLen=[0]*500
    LinkStart=[0]*500
    G=[0]*500


    for a in range (0,FwdLen):
        b=0
        LinkNumber=0

        while True:
            if(((PrimerFwd[a+b]=='A') and (PrimerRew[b]=='T'))  or ((PrimerFwd[a+b]=='T') and (PrimerRew[b]=='A'))  or ((PrimerFwd[a+b]=='G') and (PrimerRew[b]=='C'))  or ((PrimerFwd[a+b]=='C') and (PrimerRew[b]=='G'))):
                LinkOpen=True
                LinkExist=True
                LinkLen[LinkNumber]=1
                LinkStart[LinkNumber]=a+b
                b+=1
                while(((PrimerFwd[a+b]=='A') and (PrimerRew[b]=='T'))  or ((PrimerFwd[a+b]=='T') and (PrimerRew[b]=='A'))  or ((PrimerFwd[a+b]=='G') and (PrimerRew[b]=='C'))  or ((PrimerFwd[a+b]=='C') and (PrimerRew[b]=='G'))):
                    LinkLen[LinkNumber]=LinkLen[LinkNumber]+1
                    b+=1

                LinkNumber+=1

            else:
                b+=1
                if(b>=min(RewLen,FwdLen-a)):
                    break




        if(LinkNumber):
            for c in range(0,LinkNumber):
                if(LinkLen[c]>1):
                    filetext=[0]*500
                    for d in range (0,LinkLen[c]):
                        filetext[d]=PrimerFwd[LinkStart[c]+d]

                    filetext=filetext[:LinkLen[c]]

                    H=Entalpy(filetext)
                    S=Entropy(filetext)
                    G[c]=H-DimersTemperature*S

                else:
                    if((PrimerFwd[LinkStart[c]]=='A')  or (PrimerFwd[LinkStart[c]]=='T')):
                        H=-3717
                        S=-10.65
                        G[c]=(H-DimersTemperature*S)/2
                    else:
                        H=-4733
                        S=-11.91
                        G[c]=(H-DimersTemperature*S)/2

                    if(SaltSwitcher):
                        G[c]=G[c]+16.6*S*log(SaltConc,2.71828)/exp(1)



            RMid=0
            e=0
            LinkMidSize=0
            for c in range(0,LinkNumber):
                for d in range(0,LinkLen[c]):
                    RMid=RMid+float(LinkStart[c]+d)
                    e+=1
                LinkMidSize=LinkMidSize+float(LinkLen[c])

            RMid=RMid/e
            LinkMidSize=LinkMidSize/(float(LinkNumber))
            e=0
            RModMid=0
            for c in range (0,LinkNumber):
                for d in range(0,LinkLen[c]):
                    RModMid=RModMid+max(RMid,float(LinkStart[c]+d))-min(RMid,float(LinkStart[c]+d))
                    e+=1

            RModMid=RModMid/e
            if(RModMid<1):
                RModMid=1

            if(RModMid!=0):
                Expo=exp(log(RModMid/(float(max(LinkMidSize-1,1.0))),2.71828)*0.75)*1.5
            else:
                Expo=RModMid
            En=0

            for c in range(0,LinkNumber):
                if(G[c]<0):
                    En=En+exp(Expo*log((-1)*G[c],2.71828))

                else:
                    En=En-exp(Expo*log(G[c],2.71828))


            if(En>0):
                En=(-1)*exp(log(En)/Expo)
            else:
                En=exp(log((-1)*En)/Expo)


            if(a+RewLen>FwdLen):
                EnT=En*(FPen/2)*(float(RewLen-(FwdLen-a)))/(max((float(RewLen)),5.0))

            EnF=En*(FPen/2)*(float(a))/(max((float(FwdLen)),5.0))
            En=En+EnT+EnF
            if(((-1)*En)>MaxEn):
                MaxEn=(-1)*En
                if(WriteSwitcher):
                    for b in range(0,FwdLen):
                        lower[b]=PrimerFwd[b]

                    lowerLen=FwdLen
                    for b in range(0,a+1):
                        upper[b]=' '

                    for b in range(0,RewLen):
                        upper[a+b]=PrimerRew[b]

                    upperLen=a+RewLen

                    for b in range(0,FwdLen):
                        if(((upper[b]=='A') and (lower[b]=='T'))  or ((upper[b]=='T') and (lower[b]=='A'))  or ((upper[b]=='G') and (lower[b]=='C'))  or ((upper[b]=='C') and (lower[b]=='G'))):
                            binds[b]='|'
                        else:
                            binds[b]=' '


                    bindsLen=min(FwdLen,a+RewLen)

    if LinkExist:
        DimerInfo=TDimerInfo()
        DimerInfo.upper=''.join(upper[:upperLen])
        DimerInfo.lower=''.join(lower[:lowerLen])
        DimerInfo.binds=''.join(binds[:bindsLen])
        DimerInfo.Energy=-MaxEn/1000/4
        return DimerInfo
    else:
        return None

SaltSwitcher=False
ThermoMethodics=1
DimersTemperature=298
WriteSwitcher=True
SaltConc=0.05
PrimerConc=5e-8

