class TSourceDataset_Info():
    def __init__(self,Dataset_Abbreviation):
        self.Dataset_Abbreviation = Dataset_Abbreviation
        self.Methyl_samples_count_T = float('NaN')
        self.Methyl_samples_count_C = float('NaN')
        self.Methyl_samples_count_paired = float('NaN')

        self.RNASeq_samples_count_T = float('NaN')
        self.RNASeq_samples_count_C = float('NaN')
        self.RNASeq_samples_count_paired =float('NaN')

        self.RNASeq_and_Methyl_samples_count_T = float('NaN')
        self.RNASeq_and_Methyl_samples_count_C = float('NaN')
        self.RNASeq_and_Methyl_samples_count_paired = float('NaN')

        self.miRNASeq_samples_count_T = float('NaN')
        self.miRNASeq_samples_count_C = float('NaN')
        self.miRNASeq_samples_count_paired = float('NaN')

        self.RNASeq_and_miRNASeq_samples_count_T = float('NaN')
        self.RNASeq_and_miRNASeq_samples_count_C = float('NaN')
        self.RNASeq_and_miRNASeq_samples_count_paired = float('NaN')

        TCGA_Codes = {
            'LAML': 'Acute Myeloid Leukemia',
            'ACC': 'Adrenocortical carcinoma',
            'BLCA': 'Bladder Urothelial Carcinoma',
            'LGG': 'Brain Lower Grade Glioma',
            'BRCA': 'Breast invasive carcinoma',
            'CESC': 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
            'CHOL': 'Cholangiocarcinoma',
            'COAD': 'Colon adenocarcinoma',
            'ESCA': 'Esophageal carcinoma',
            'FPPP': 'FFPE Pilot Phase II',
            'GBM': 'Glioblastoma multiforme',
            'HNSC': 'Head and Neck squamous cell carcinoma',
            'KICH': 'Kidney Chromophobe',
            'KIRC': 'Kidney renal clear cell carcinoma',
            'KIRP': 'Kidney renal papillary cell carcinoma',
            'LIHC': 'Liver hepatocellular carcinoma',
            'LUAD': 'Lung adenocarcinoma',
            'LUSC': 'Lung squamous cell carcinoma',
            'DLBC': 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
            'MESO': 'Mesothelioma',
            'OV': 'Ovarian serous cystadenocarcinoma',
            'PAAD': 'Pancreatic adenocarcinoma',
            'PCPG': 'Pheochromocytoma and Paraganglioma',
            'PRAD': 'Prostate adenocarcinoma',
            'READ': 'Rectum adenocarcinoma',
            'SARC': 'Sarcoma',
            'SKCM': 'Skin Cutaneous Melanoma',
            'STAD': 'Stomach adenocarcinoma',
            'TGCT': 'Testicular Germ Cell Tumors',
            'THYM': 'Thymoma',
            'THCA': 'Thyroid carcinoma',
            'UCS': 'Uterine Carcinosarcoma',
            'UCEC': 'Uterine Corpus Endometrial Carcinoma',
            'UVM': 'Uveal Melanoma'
        }

        self.Dataset_Name = None
        if self.Dataset_Abbreviation in TCGA_Codes:
            self.Dataset_Name = TCGA_Codes[self.Dataset_Abbreviation]




class TRNASeq_BasicInfo():
    def __init__(self,GeneName):
        self.GeneName = GeneName
        self.LogFC_paired = float('NaN')
        self.FDR_paired = float('NaN')
        self.LogFC_pooled = float('NaN')
        self.FDR_pooled = float('NaN')

class TCpG_BasicInfo():
    def __init__(self):
        # self.CpG_id = CpG_id
        # self.Chr = ''
        # self.Pos = float('NaN')
        self.HyperMeth_score = float('NaN')
        self.HypoMeth_score = float('NaN')
        self.Rs_paired = dict()
        self.Rs_pooled = dict()