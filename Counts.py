import pandas as pd
import numpy as np

class Domain(object):
    """A class for each domain: Archaea, Bacteria, Eukaryotes, Viruses, Unknown
    """

    domDict = {}

    def __init__(self, name):
        """ Initializer for a Domain"""

        # print('Initializing a new Domain')
        self.name = name
        self.species = {}
        self.isoforms = {}
        
        

class Species(object):
    """
    A class for each species, containing domain information, isoform information, and amino acid information
    """
    speciesDict = {}
    def __init__(self, name):
        """
        Initializer for a Species
        """

        self.name = name
        self.isoforms = {}
        self.totalAA = 0
        Species.speciesDict[name] = self
        aaSet = ['ala','arg', 'asn','asp','cys','glu', 'gln', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']
        for aa in aaSet:
            setattr(self, aa, 0)
    

class Protein(object):
    """
    A class for each protein, containing its data, id, protein#, species, domain, and amino acid information
    """
    proteinDict = {}
    proteinAACounts = []
    proteinAACountsArray = None
    AllProteinsAA = 0
    def __init__(self, id, p, asterisk, aaSeq):
        """
        Initializer for a new protein
        """

        self.id = id
        self.p = p
        self.asterisk = asterisk
        self.aaSeq = aaSeq
        self.domain = None
        self.species = None
       
        self.aaCount()
        if id in Protein.proteinDict:
            print("\n\n\nDUPLICATE PROTEINNNNN", id)
        Protein.proteinDict[id] = self
        

    def aaCount(self):
        """
        Creates a class variable for each amino abbreviation, then counts 
        the number of amino acids in a protein's amino acid sequence and saves the total to the
        respective variable, such as Protein.arg = 12
        Also saves the length of the amino acid sequence, such as Protein.len = 200
        """
        #print('Beginning aaCount for protein ['+self.id+']')
        len = 0
        totalAA = 0

        ala = 0 
        arg = 0
        asn = 0
        asp = 0
        cys = 0
        glu = 0
        gln = 0
        gly = 0
        his = 0
        ile = 0
        leu = 0
        lys = 0
        met = 0
        phe = 0
        pro = 0
        ser = 0
        thr = 0
        trp = 0
        tyr = 0
        val = 0

        for i in self.aaSeq:
            len += 1
            if i == 'A':
                ala += 1
                totalAA += 1
            elif i == 'R':
                arg += 1
                totalAA += 1
            elif i == 'N':
                asn += 1
                totalAA += 1
            elif i == 'D':
                asp += 1
                totalAA += 1
            elif i == 'C':
                cys += 1
                totalAA += 1
            elif i == 'E':
                glu += 1
                totalAA += 1
            elif i == 'Q':
                gln += 1
                totalAA += 1
            elif i == 'G':
                gly += 1
                totalAA += 1
            elif i == 'H':
                his += 1
                totalAA += 1
            elif i == 'I':
                ile += 1
                totalAA += 1
            elif i == 'L':
                leu += 1
                totalAA += 1
            elif i == 'K':
                lys += 1
                totalAA += 1
            elif i == 'M':
                met += 1
                totalAA += 1
            elif i == 'F':
                phe += 1
                totalAA += 1
            elif i == 'P':
                pro += 1
                totalAA += 1
            elif i == 'S':
                ser += 1
                totalAA += 1
            elif i == 'T':
                thr += 1
                totalAA += 1
            elif i == 'W':
                trp += 1
                totalAA += 1
            elif i == 'Y':
                tyr += 1
                totalAA += 1
            elif i == 'V':
                val += 1
                totalAA += 1

        if totalAA > 0:
            alaPercent = ala  / totalAA * 100
            argPercent = arg  / totalAA * 100
            asnPercent = asn  / totalAA * 100
            aspPercent = asp  / totalAA * 100
            cysPercent = cys  / totalAA * 100
            gluPercent = glu  / totalAA * 100
            glnPercent = gln  / totalAA * 100
            glyPercent = gly  / totalAA * 100
            hisPercent = his  / totalAA * 100
            ilePercent = ile  / totalAA * 100
            leuPercent = leu  / totalAA * 100
            lysPercent = lys  / totalAA * 100
            metPercent = met  / totalAA * 100
            phePercent = phe  / totalAA * 100
            proPercent = pro  / totalAA * 100
            serPercent = ser  / totalAA * 100
            thrPercent = thr  / totalAA * 100
            trpPercent = trp  / totalAA * 100
            tyrPercent = tyr  / totalAA * 100
            valPercent = val  / totalAA * 100
        else: 
            alaPercent = 0
            argPercent = 0
            asnPercent = 0
            aspPercent = 0
            cysPercent = 0
            gluPercent = 0
            glnPercent = 0
            glyPercent = 0
            hisPercent = 0
            ilePercent = 0
            leuPercent = 0
            lysPercent = 0
            metPercent = 0
            phePercent = 0
            proPercent = 0
            serPercent = 0
            thrPercent = 0
            trpPercent = 0
            tyrPercent = 0
            valPercent = 0

        data = (self.id, len, totalAA, ala,arg, asn,asp,cys,glu, gln, gly, his,ile,leu,lys,met,phe,pro, ser, thr, trp, tyr,val, alaPercent, argPercent, asnPercent, aspPercent, cysPercent, gluPercent, glnPercent, glyPercent, hisPercent, ilePercent, leuPercent, lysPercent, metPercent, phePercent, proPercent, serPercent, thrPercent, trpPercent, tyrPercent, valPercent)
        Protein.proteinAACounts.append(data)

        self.len = 0
        self.totalAA = 0

        self.ala = 0 
        self.arg = 0
        self.asn = 0
        self.asp = 0
        self.cys = 0
        self.glu = 0
        self.gln = 0
        self.gly = 0
        self.his = 0
        self.ile = 0
        self.leu = 0
        self.lys = 0
        self.met = 0
        self.phe = 0
        self.pro = 0
        self.ser = 0
        self.thr = 0
        self.trp = 0
        self.tyr = 0
        self.val = 0
        for i in self.aaSeq:
            self.len += 1
            if i == 'A':
                self.ala += 1
                self.totalAA += 1
            elif i == 'R':
                self.arg += 1
                self.totalAA += 1
            elif i == 'N':
                self.asn += 1
                self.totalAA += 1
            elif i == 'D':
                self.asp += 1
                self.totalAA += 1
            elif i == 'C':
                self.cys += 1
                self.totalAA += 1
            elif i == 'E':
                self.glu += 1
                self.totalAA += 1
            elif i == 'Q':
                self.gln += 1
                self.totalAA += 1
            elif i == 'G':
                self.gly += 1
                self.totalAA += 1
            elif i == 'H':
                self.his += 1
                self.totalAA += 1
            elif i == 'I':
                self.ile += 1
                self.totalAA += 1
            elif i == 'L':
                self.leu += 1
                self.totalAA += 1
            elif i == 'K':
                self.lys += 1
                self.totalAA += 1
            elif i == 'M':
                self.met += 1
                self.totalAA += 1
            elif i == 'F':
                self.phe += 1
                self.totalAA += 1
            elif i == 'P':
                self.pro += 1
                self.totalAA += 1
            elif i == 'S':
                self.ser += 1
                self.totalAA += 1
            elif i == 'T':
                self.thr += 1
                self.totalAA += 1
            elif i == 'W':
                self.trp += 1
                self.totalAA += 1
            elif i == 'Y':
                self.tyr += 1
                self.totalAA += 1
            elif i == 'V':
                self.val += 1
                self.totalAA += 1

        if self.totalAA > 0:
            self.alaPercent = self.ala  / self.totalAA * 100
            self.argPercent = self.arg  / self.totalAA * 100
            self.asnPercent = self.asn  / self.totalAA * 100
            self.aspPercent = self.asp  / self.totalAA * 100
            self.cysPercent = self.cys  / self.totalAA * 100
            self.gluPercent = self.glu  / self.totalAA * 100
            self.glnPercent = self.gln  / self.totalAA * 100
            self.glyPercent = self.gly  / self.totalAA * 100
            self.hisPercent = self.his  / self.totalAA * 100
            self.ilePercent = self.ile  / self.totalAA * 100
            self.leuPercent = self.leu  / self.totalAA * 100
            self.lysPercent = self.lys  / self.totalAA * 100
            self.metPercent = self.met  / self.totalAA * 100
            self.phePercent = self.phe  / self.totalAA * 100
            self.proPercent = self.pro  / self.totalAA * 100
            self.serPercent = self.ser  / self.totalAA * 100
            self.thrPercent = self.thr  / self.totalAA * 100
            self.trpPercent = self.trp  / self.totalAA * 100
            self.tyrPercent = self.tyr  / self.totalAA * 100
            self.valPercent = self.val  / self.totalAA * 100
        else: 
            self.alaPercent = 0
            self.argPercent = 0
            self.asnPercent = 0
            self.aspPercent = 0
            self.cysPercent = 0
            self.gluPercent = 0
            self.glnPercent = 0
            self.glyPercent = 0
            self.hisPercent = 0
            self.ilePercent = 0
            self.leuPercent = 0
            self.lysPercent = 0
            self.metPercent = 0
            self.phePercent = 0
            self.proPercent = 0
            self.serPercent = 0
            self.thrPercent = 0
            self.trpPercent = 0
            self.tyrPercent = 0
            self.valPercent = 0

        Protein.AllProteinsAA += self.totalAA

    def clearProteinAACounts():
        """
        Clears Protein.proteinAACounts
        """
        Protein.proteinAACounts = []

class Sample(object):
    """
    A class for each sample, containing its name and amino acid sums
    """
    sampleDict = {}
    abundanceData = None

    def __init__(self, name):
        """
        Initializes a Sample object, assigns its name.
        name is a string
        proteinAbundance is a Series with protein ids as the index and floats as the content.
        aaData is a DataFrame containing total aa counts and percentages for each domain.
        """
        self.name = name
        self.proteinAbundance = None
        self.aaData = pd.DataFrame(np.full(shape = (20,12), fill_value = 0),index = ['ala','arg', 'asn','asp','cys','glu', 'gln', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val'], columns = ['Total', 'Percentage', 'Total_Archaea', 'Total_Bacteria', 'Total_Eukaryota', 'Total_Viruses', 'Total_Unknown','Percent_Archaea','Percent_Bacteria','Percent_Eukaryota','Percent_Viruses','Percent_Unknown'],dtype= 'float64')
        self.aaData.name = name

        