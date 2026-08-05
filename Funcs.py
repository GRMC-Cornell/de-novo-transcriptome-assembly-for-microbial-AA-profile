from Counts import *
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 10000)
pd.set_option('display.max_colwidth', None)

import matplotlib.pyplot as plt
import os


def isRegList(lst):
    """
    Returns True if lst is a two-dimensional non-ragged list of strings, otherwise returns False.

    lst can be anything.
    """
    #print('Beginning isRegList')
    lengths = []
    if isinstance(lst, list):
        for row in lst:
            if isinstance(row, list):
                lengths = lengths + [len(row)]
                for col in row:
                    if isinstance(col, str):
                        pass
                    else:
                        return False
            else:
                return False
        width = len(lst[0])
        if min(lengths) == width and max(lengths) == width:
            return True
        else:
            return False
    else:
        return False


def isCatalogue(dct):
    """
    Returns True if dct is a dict in the form {'string', ('string','string'), ...}, where the strings are not empty, 
    otherwise returns False.
    """
    #print('Beginning isCatalogue')
    if not isinstance(dct, dict):
        return False
    else:
        if not len(dct) > 0:
            return False
        else:
            for key, value in dct.items():
                if not isinstance(key, str):
                    return False
                elif not isinstance(value, tuple):
                    return False
                elif not len(value) == 2:
                    return False
                else:
                    for i in value:
                        if not isinstance(i, str):
                            return False
                        if len(i) == 0:
                            return False
    return True


def isProtEntry(string):
    """
    Returns True if string has the form 'ProteinID.isoform ...)aminoacidsequence*' where ProteinID can be any 
    combination of letters (either case), numbers, and underscores starting with a letter or underscore, isoform must be the form 'p#' where # is any 
    combination of digits, ... can be anything that starts with a ' ' space, and aminoacidsequence must be only letters. 
    The ending asterisk is optional.
    Otherwise, returns False.

    string can be anything.
    """
    #print('Beginning isProtEntry')
    if isinstance(string, str):
        if string.find('.') > 0:
            if string[:string.find('.')].isidentifier():
                if string[string.find('.')+1] == 'p':
                    if string[string.find('.')+2:string.find(' ')].isdecimal():
                        if string.rfind(')') < len(string) and string.rfind(')') > string.find(' '):
                            if string[-1] == '*':
                                if string[string.rfind(')')+1:-1].isalpha():
                                    return True
                            else:
                                if string[string.rfind(')')+1:].isalpha():
                                    return True
    return False


def isProtMatch(data):
    """
    Returns True if data is a dict with the form {ProteinID: ProteinObject, ...} where ProteinID 
    is a string and ProteinObject is a Protein Object.
    Returns False otherwise

    data can be anything
    """
    if isinstance(data, dict):
        if len(data) > 0:
            for key, value in data.items():
                if isinstance(key, str):
                    if isinstance(value, Protein):
                        pass
                    else:
                        return False
                else:
                    return False
            return True
    return False


def createDomains():
    """
    Initializes the five domain classes - Archaea, Bacteria, Eukaryota, Viruses, Unknown 
    into global variables of their respective names.
    """
   #print('Beginning createDomains')

    Domain.domDict['Archaea'] = Domain('Archaea')
    Domain.domDict['Bacteria'] = Domain('Bacteria')
    Domain.domDict['Eukaryota'] = Domain('Eukaryota')
    Domain.domDict['Viruses'] = Domain('Viruses')
    Domain.domDict['Unknown'] = Domain('Unknown')
    # print('Created Domains')
    

def clearDomains():
    """
    Removes all domain information from memory, use for testing.
    """
    #print('Beginning clearDomains')
    Domain.domDict = {}
    # print('Cleared Domains')


def extractBLAST(filename):
    """
    Returns data from Taxonomy file as a pandas DataFrame.

    filename is a .csv file with entries separated by newlines, with commas separating values.
    """
   #print('Beginning extractBLAST')
    assert type(filename) == str, 'filename is not a string'
    assert len(filename) > 0, 'filename is empty'
    assert filename[-4:] == '.csv', 'filename is not a .csv file'

    df = pd.read_csv(filename)
    #print('BLASTlist is: '+ str(BLASTlist))    
    return df


def taxParse(data):
    """
    Returns a dict of all proteins assigned to a list of their Domain and Species classification.
    in the form {id: (domain, species), ...} where domain is 'Unknown' and species is 'None' if the information does not exist.

    data is a pandas DataFrame with entries as a string or nan, where the first column contains the proteinIDs, and there is a column named 'Domain' 
    and a column named 'Species'.
    """
   #print('Beginning taxParse')
    assert isinstance(data, pd.DataFrame), 'data is not a pandas DataFrame'
    assert len(data) > 0, 'DataFrame is empty'
    assert 'Domain' in data.columns, "data does not have a column named 'Domain'"
    assert 'Species' in data.columns, "data does not have a column named 'Species'"

    catalogue = {}
    data = data.fillna('None')
    domList = list(data['Domain'])
    assert all(isinstance(i, str) for i in domList)
    specList = list(data['Species'])
    assert all(isinstance(i, str) for i in specList)
    idList = list(data[data.columns[0]])
    assert all(isinstance(i, str) for i in idList)
    tax = list(zip(domList,specList))
    catalogue = dict(zip(idList,tax))

    return catalogue


def assignSpecies(dct):
    """
    Creates and assigns Species objects to the domains.
    Takes a dictionary of proteins in the form {'Protein Name': ['domain', 'species'],...} where domain can be 'Unknown'and species can be 'None'.
    If species is 'None', then no new species object will be created
    
    dct is a dict of protein keys with a value in the form {'Protein Name': ['domain', 'species'],...} where Protein Name, domain, and species
    are all strings of any length greater than 0. The value must have a length of 2
    """
   #print('Beginning assignSpecies')
    assert isCatalogue(dct)

    for key, value in dct.items():
        if value[0] == 'Archaea':
            if value[1] != 'None':
                if value[1] not in Domain.domDict['Archaea'].species:
                    Domain.domDict['Archaea'].species[value[1]] = Species(value[1])
                    # print('Added '+value[1]+' to Archaea species list')
        elif value[0] == 'Bacteria':
            if value[1] != 'None':
                if value[1] not in Domain.domDict['Bacteria'].species:
                    Domain.domDict['Bacteria'].species[value[1]] = Species(value[1])
                    # print('Added '+value[1]+' to Bacteria species list')
        elif value[0] == 'Eukaryota':
            if value[1] != 'None':
                if value[1] not in Domain.domDict['Eukaryota'].species:
                    Domain.domDict['Eukaryota'].species[value[1]] = Species(value[1])
                    # print('Added '+value[1]+' to Eukaryota species list')
        elif value[0] == 'Viruses':
            if value[1] != 'None':
                if value[1] not in Domain.domDict['Viruses'].species:
                    Domain.domDict['Viruses'].species[value[1]] = Species(value[1])
                    # print('Added '+value[1]+' to Viruses species list')
        elif value[0] == 'Unknown':
            if value[1] != 'None':
                if value[1] not in Domain.domDict['Unknown'].species:
                    Domain.domDict['Unknown'].species[value[1]] = Species(value[1])
                    # print('Added '+value[1]+' to Unknown species list')


def extractDataset(filename):
    """
    Returns data from filename as a single string, without whitespace.

    filename is a string of a .pep file.
    """
   #print('Beginning extractDataset')
    assert isinstance(filename, str), "filename is not a string"
    assert len(filename) > 0, "filename is empty"
    assert filename[-4:] == '.pep', "file type is not '.pep'"

    try:
        with open(filename, "r") as dataset:
            datasetContent = ""
            for n in dataset:
                datasetContent = datasetContent + n.strip()
        return datasetContent
    except FileNotFoundError:
        raise FileNotFoundError("filename does not exist.")


def transformDataset(dataset):
    """
    Returns a dictionary of protein ids matched with their objects in the form {ProteinID: Object, ...}.
    ProteinID is Isoform.Protein#
    Parses a string of protein entries separated by '>' into individual strings, then parses individual strings for 
    relevant information, creating a Protein object for each and adding it to a dictionary.

    dataset is a string of protein information with correct format, with each protein entry beginning with '>'.
    """
    ### NEEDS OPTIMIZATION ###
   #print('Beginning transformDataset')

    assert isinstance(dataset, str), "dataset is not a string"
    assert len(dataset) > 0, "dataset is empty"
    assert dataset[0] == '>', "First protein entry does not start with '>'."

    entries = dataset.split('>')[1:]
    protMatch = {}
    for protein in entries:
        assert isProtEntry(protein)
        new = newProtein(protein)
        protMatch[new[0]] = new[1]
        #print('protMatch is now '+str(protMatch))
    # All protein objects have been made, turn Protein.proteinAACounts into an array
    Protein.proteinAACountsArray = np.array(Protein.proteinAACounts, dtype = np.dtype([('id','<U50'),('length',np.float64),('totalAA', np.float64),('ala', np.float64),('arg', np.float64),('asn', np.float64),('asp', np.float64),('cys', np.float64),('gln', np.float64),('glu', np.float64),('gly', np.float64),('his', np.float64),('ile', np.float64),('leu', np.float64),('lys', np.float64),('met', np.float64),('phe', np.float64),('pro', np.float64),('ser', np.float64),('thr', np.float64),('trp', np.float64),('tyr', np.float64),('val', np.float64),('alaPercent', np.float64),('argPercent', np.float64),('asnPercent', np.float64),('aspPercent', np.float64),('cysPercent', np.float64),('glnPercent', np.float64),('gluPercent', np.float64),('glyPercent', np.float64),('hisPercent', np.float64),('ilePercent', np.float64),('leuPercent', np.float64),('lysPercent', np.float64),('metPercent', np.float64),('phePercent', np.float64),('proPercent', np.float64),('serPercent', np.float64),('thrPercent', np.float64),('trpPercent', np.float64),('tyrPercent', np.float64),('valPercent', np.float64)]))
    # print('ProteinAACounts Array is \n')
    # print(Protein.proteinAACountsArray)
    assert len(entries) == len(protMatch), "length of entries and protMatch are not equal"
    
    return protMatch


def newProtein(data):
    """
    Helper for transformDataset()
    Creates a new protein object from data, parses isoform, protein#, asterisk status, and amino acid sequence.
    Returns the new protein matched with its name, as a list: ['proteinID', Protein-Object] where proteinID is Isoform.protein#

    data is a non-empty string containing a single protein entry with id, isoform, and amino acid sequence.
    data has the form 'ProteinID.isoform ...)aminoacidsequence*' where ProteinID can be any combination of letters (either case), 
    numbers, and underscores starting with a letter or underscore, isoform must be the form 'p#' where # is any combination of 
    digits, ... can be anything that starts with a ' ' space, and aminoacidsequence must be only letters. The ending asterisk is optional.
    """
    #print('Beginning newProtein')

    assert isinstance(data, str), "data is not a string"
    assert len(data) > 0, "data is empty"
    assert isProtEntry(data), "data has an invalid format"

    id = data[:data.find(' ')]
    p = id[id.find('.p')+2:]
    if data.endswith('*'):
        asterisk = True
        aaSeq = data[data.rfind(')')+1:-1]
    else:
        asterisk = False
        aaSeq = data[data.rfind(')')+1:]
    # print('ID is: '+id)
    # print('protein is: '+p)
    # print('Asterisk is: '+str(asterisk))
    # print('aaSeq is: '+aaSeq)
    return [id, Protein(id, p, asterisk, aaSeq)]
        

def assignIsoform(protMatch, catalogue):
    """
    Assigns protein objects to their respective species objects and domain objects.

    protMatch is a dict of proteins in the form {id : object}
    catalogue is a dict of proteins in the form {id : [domain, species]}
    Every key in protMatch must have a corresponding key in catalogue.
    """
   #print('Beginning assignIsoform')
    assert isProtMatch(protMatch), "protMatch as invalid format"
    assert isCatalogue(catalogue), "catalogue has invalid format"


    for key, value in protMatch.items():
        id = key
        isoform = key[: key.index('.')]
        #print('isoform is: '+ isoform)
        if isoform not in catalogue:
            Domain.domDict['Unknown'].isoforms[key] = value
            value.domain = 'Unknown'
            value.species = 'None'
            #print('isoform '+isoform+ ' not in catalogue, added to Unknown domain')
        else:
            domain = catalogue[isoform][0]
            #print('domain is '+domain)
            species = catalogue[isoform][1]
            #print('species is '+ species)
            if domain != 'None':
                value.domain = domain
            elif domain =='None':
                value.domain = 'Unknown'
            value.species = species
            if domain == 'Archaea':
                Domain.domDict['Archaea'].isoforms[id] = value
                #print('Added '+id+' to Archaea')
                if species != 'None':
                    Domain.domDict['Archaea'].species[species].isoforms[id] = value
                    #print('Added '+id+' to '+species+' in Archaea')
            elif domain == 'Bacteria':
                Domain.domDict['Bacteria'].isoforms[id] = value
                #print('Added '+id+' to Bacteria')
                if species != 'None':
                    Domain.domDict['Bacteria'].species[species].isoforms[id] = value
                    #print('Added '+id+' to '+species+' in Bacteria')
            elif domain == 'Eukaryota':
                Domain.domDict['Eukaryota'].isoforms[id] = value
                #print('Added '+id+' to Eukaryota')
                if species != 'None':
                    Domain.domDict['Eukaryota'].species[species].isoforms[id] = value
                    #print('Added '+id+' to '+species+' in Eukaryota')
            elif domain == 'Viruses':
                Domain.domDict['Viruses'].isoforms[id] = value
                #print('Added '+id+' to Viruses')
                if species != 'None':
                    Domain.domDict['Viruses'].species[species].isoforms[id] = value
                    #print('Added '+id+' to '+species+' in Viruses')
            elif domain == 'None':
                Domain.domDict['Unknown'].isoforms[id] = value
                #print('Added '+id+' to Unknown')
                if species != 'None':
                    Domain.domDict['Unknown'].species[species].isoforms[id] = value
                    #print('Added '+id+' to '+species+' in Unknown')


def extractAbundance(filename):
    """
    Opens filename using pandas and creates a DataFrame. Assigns DataFrame to Sample.abundanceData.

    filename is a nonempty string of a valid .csv file where the first line contains headers to the columns.
    """
   #print('Beginning extractAbundance()')
    assert isinstance(filename, str), 'filename is not a string'
    assert len(filename) > 0, 'filename is empty'
    assert filename[-4:] == '.csv', 'filename is not a .csv file'
    try:
        abundanceData = pd.read_csv(filename, sep='	', index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError("No such file exists")
    abundanceData.columns = abundanceData.columns.str.strip()
    Sample.abundanceData = abundanceData
    return abundanceData


def createSamples(abundanceData):
    """
    Creates a Sample object for every sample in the abundanceData DataFrame

    abundanceData is a pandas DataFrame with columns as sample names. Each column header must be a non-empty string.
    """
   #print('Beginning createSamples()')
    assert isinstance(abundanceData, pd.DataFrame), "abundanceData is not a DataFrame"
    assert len(abundanceData) > 0, "abundanceData is empty"

    for s in abundanceData.columns:
        assert isinstance(s, str), 'Column '+str(s)+' is not a string'
        assert len(s) > 0, 'Empty column'

        Sample.sampleDict[s] = Sample(s)
    Sample.sampleDict = dict(sorted(Sample.sampleDict.items()))
    #print('sampleDict is '+ str(Sample.sampleDict))
    

def applySampleMetadata():
    """
    Categorizes each sample as feed inefficient vs feed efficient based on metadata.
    """
    for sampObj in Sample.sampleDict.values():
        if sampObj.name == '18048XR-81-01_S0_L001':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-02_S1_L004':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-03_S2_L003':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-04_S3_L003':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-05_S4_L003':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-06_S5_L003':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-07_S6_L003':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-08_S7_L003':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-09_S0_L001':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-10_S3_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-11_S8_L003':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-12_S0_L001':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-13_S5_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-14_S0_L001':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-15_S0_L001':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-16_S6_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-17_S7_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-18_S8_L004':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-19_S9_L004':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-20_S10_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-21_S11_L004':
            sampObj.feedEfficient = True
        elif sampObj.name == '18048XR-81-22_S12_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-23_S13_L004':
            sampObj.feedEfficient = False
        elif sampObj.name == '18048XR-81-24_S11_L003':
            sampObj.feedEfficient = True
        

def clearSamples():
    """
    Clears all samples objects from the class dictionary.
    """
    Sample.sampleDict = {}


def differentiateAbundance(abundanceData):
    """
    Breaks up abundanceData DataFrame into individual pandas Series by sample (by column).
    Assigns each Series to its respective Sample object as Sample.proteinAbundance

    abundanceData is a non-empty pandas DataFrame where the columns are sample names that exist in sampleDict.
    """
   #print('Beginning differentiateAbundance()')
    assert isinstance(abundanceData,pd.DataFrame)
    assert len(abundanceData) > 0

    for col in abundanceData:
        assert isinstance(col, str)
        assert col in Sample.sampleDict
        if col in Sample.sampleDict:

            Sample.sampleDict[col].proteinAbundance = abundanceData[col]
        

def computeRawAA(protMatch):
    """
    Creates a pandas DataFrame for each Sample in Sample.sampleDict, containing total AA count information. 
    Stores it in the attribute Sample.aaData
    Before calling, Sample.sampleDict must not be empty.

    protMatch is a non-empty dictionary of protein IDs (strings) matched with Protein objects.

    sample is a string of a sample name (column in the abundance file)
    """
    ### NEEDS OPTIMIZATION ###
   #print('beginning computeRawAA()')
    assert isProtMatch(protMatch), 'protMatch as invalid format'
    assert len(Sample.sampleDict) > 0, 'Sample.sampleDict is empty'

    aaSet = ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']
    directory = 'Samples_IsoformData'
    os.makedirs(directory, exist_ok=True)
    for sampObj in Sample.sampleDict.values():
        full_path = os.path.join(directory, f'{sampObj.name}_IsoformData.csv')
        with open(full_path, 'w') as file:
            for aa in aaSet:
                file.write(f'{aa},')
            file.write('Transcript\n')
    tempSampData = {}
    abundanceData = Sample.abundanceData.to_dict(orient = 'index')
    #print('abundanceData is:\n'+str(abundanceData))
    for sample in Sample.sampleDict.keys():
        tempSampData[sample] = [[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]]

    totalAA = 0
    for protObj in protMatch.values():
        #print('starting new protein\n')
        isoform = protObj.id[:protObj.id.find('.')]
        #print('isoform is: '+isoform)
        protAbundance = abundanceData[isoform]
        #print('protAbundance is:\n'+str(protAbundance))
        aaValues = []
        for aa in aaSet:
            aaValues.append(getattr(protObj, aa))
        for sampObj in Sample.sampleDict.values():
            #print('starting new sample for protein')
            mult = protAbundance[sampObj.name]
            #print('mult is '+str(mult))
            full_path = os.path.join(directory, f'{sampObj.name}_IsoformData.csv')
            with open(full_path, 'a') as file:
                for n in range(20):
                    aaVal = aaValues[n]
                    file.write(str(aaVal*mult)+',')
                    tempSampData[sampObj.name][n][0] += aaVal * mult
                    if protObj.domain == 'Archaea':
                        tempSampData[sampObj.name][n][2] += aaVal * mult
                    elif protObj.domain == 'Bacteria':
                        tempSampData[sampObj.name][n][3] += aaVal * mult
                    elif protObj.domain == 'Eukaryota':
                        tempSampData[sampObj.name][n][4] += aaVal * mult
                    elif protObj.domain == 'Viruses':
                        tempSampData[sampObj.name][n][5] += aaVal * mult
                    elif protObj.domain == 'Unknown':
                        tempSampData[sampObj.name][n][6] += aaVal * mult
                file.write(f'{protObj.id}\n')
    
    Domain.totalAA = totalAA                  
    for sampObj in Sample.sampleDict.values():
        sampObj.aaData = pd.DataFrame(tempSampData[sampObj.name], columns = ['Total',  'Percentage',  'Total_Archaea',  'Total_Bacteria',  'Total_Eukaryota',  'Total_Viruses',  'Total_Unknown',  'Percent_Archaea',  'Percent_Bacteria',  'Percent_Eukaryota',  'Percent_Viruses',  'Percent_Unknown'], dtype = np.float64, index = aaSet).round()


def computePercentAA():
    """
    Edits the Sample.aaData DataFrame within each Sample object to contain correct percentage values.
    """
   #print('Beginning computePercentAA')
    for value in Sample.sampleDict.values():
        sum = value.aaData['Total'].sum()
        for row in value.aaData.index:
            rowTotal = value.aaData.at[row, 'Total']
            #print('rowTotal is '+ str(rowTotal))
            if rowTotal != 0:
                value.aaData.at[row,'Percentage'] = round(value.aaData.at[row,'Total'] / sum * 100, 2)
                value.aaData.at[row,'Percent_Archaea'] = round(value.aaData.at[row, 'Total_Archaea'] / rowTotal * 100, 2)
                value.aaData.at[row,'Percent_Bacteria'] = round(value.aaData.at[row, 'Total_Bacteria'] / rowTotal * 100, 2)
                value.aaData.at[row,'Percent_Eukaryota'] = round(value.aaData.at[row, 'Total_Eukaryota'] / rowTotal * 100, 2)
                value.aaData.at[row,'Percent_Viruses'] = round(value.aaData.at[row, 'Total_Viruses'] / rowTotal * 100, 2)
                value.aaData.at[row,'Percent_Unknown'] = round(value.aaData.at[row, 'Total_Unknown'] / rowTotal * 100, 2)
    
    #print(Sample.sampleDict['18048XR-81-15_S0_L001'].aaData)


def computeSpeciesAACounts():
    """
    Computes a Species.totalAA value for each species object, based on the Protein.totalAA values in the Species.isoforms.
    """
   #print('Beginning computeSpeciesAACounts()')
    aaSet = ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']
    
    for obj in Species.speciesDict.values():
        tempAAData = {'ala':0,'arg':0, 'asn':0,'asp':0,'cys':0,'gln':0, 'glu':0, 'gly':0, 'his':0,'ile':0,'leu':0,'lys':0,'met':0,'phe':0,'pro':0, 'ser':0, 'thr':0, 'trp':0, 'tyr':0,'val':0}
        tempTotalAA = 0
        for protein in obj.isoforms.values():
            for aa in aaSet:
                aaVal = getattr(protein, aa)
                tempAAData[aa] = tempAAData[aa] + aaVal
                tempTotalAA += aaVal
        
        for aa in aaSet:
            setattr(obj, aa, tempAAData[aa])
            if tempTotalAA > 0:
                setattr(obj,f'{aa}Percent', round(tempAAData[aa] / tempTotalAA * 100, 2))
            else:
                setattr(obj, f'{aa}Percent', 0)
        obj.totalAA = tempTotalAA


def createExportFolder():
    """
    Creates a new folder named name to house all subsequently created files.

    name is a nonempty string. By default, name = 'Exports'.
    """
    global exportFolderName
    exportFolderName = 'GRMC_Exports'
    os.makedirs(exportFolderName, exist_ok=True)


def exportSampleData():
    """
    Creates a new file for each domain containing raw aa counts and percentage counts for each sample, with aa names as columns and sample names as rows.
    """
   #print('Begninning exportSampleData()')
    for domName, domObj in Domain.domDict.items():
        full_path = os.path.join(exportFolderName, f'{domName}_Raw_AA_Data_by_Sample.csv')
        with open(full_path, 'w') as file:
            totalAA = 0
            file.write('Sample ID, Total AA, Total ala, Total arg, Total asn, Total asp, Total cys, Total gln, Total glu, Total gly, Total his, Total ile, Total leu, Total lys, Total met, Total phe, Total pro, Total ser, Total thr, Total trp, Total tyr, Total val, Percent ala, Percent arg, Percent asn, Percent asp, Percent cys, Percent gln, Percent glu, Percent gly, Percent his, Percent ile, Percent leu, Percent lys, Percent met, Percent phe, Percent pro, Percent ser, Percent thr, Percent trp, Percent tyr, Percent val\n')
            for sampleName, sampleObj in Sample.sampleDict.items():
                totalAA += sampleObj.aaData[f'Total_{domName}'].sum()
                file.write(f"{sampleName}, {sampleObj.aaData[f'Total_{domName}'].sum()}, {sampleObj.aaData.at['ala',f'Total_{domName}']}, {sampleObj.aaData.at['arg',f'Total_{domName}']}, {sampleObj.aaData.at['asn',f'Total_{domName}']}, {sampleObj.aaData.at['asp',f'Total_{domName}']}, {sampleObj.aaData.at['cys',f'Total_{domName}']}, {sampleObj.aaData.at['glu',f'Total_{domName}']}, {sampleObj.aaData.at['gln',f'Total_{domName}']}, {sampleObj.aaData.at['gly',f'Total_{domName}']}, {sampleObj.aaData.at['his',f'Total_{domName}']}, {sampleObj.aaData.at['ile',f'Total_{domName}']}, {sampleObj.aaData.at['leu',f'Total_{domName}']}, {sampleObj.aaData.at['lys',f'Total_{domName}']}, {sampleObj.aaData.at['met',f'Total_{domName}']}, {sampleObj.aaData.at['phe',f'Total_{domName}']}, {sampleObj.aaData.at['pro',f'Total_{domName}']}, {sampleObj.aaData.at['ser',f'Total_{domName}']}, {sampleObj.aaData.at['thr',f'Total_{domName}']}, {sampleObj.aaData.at['trp',f'Total_{domName}']}, {sampleObj.aaData.at['tyr',f'Total_{domName}']}, {sampleObj.aaData.at['val',f'Total_{domName}']}, {sampleObj.aaData.at['ala','Percentage']}, {sampleObj.aaData.at['arg','Percentage']}, {sampleObj.aaData.at['asn','Percentage']}, {sampleObj.aaData.at['asp','Percentage']}, {sampleObj.aaData.at['cys','Percentage']}, {sampleObj.aaData.at['gln','Percentage']}, {sampleObj.aaData.at['glu','Percentage']}, {sampleObj.aaData.at['gly','Percentage']}, {sampleObj.aaData.at['his','Percentage']}, {sampleObj.aaData.at['ile','Percentage']}, {sampleObj.aaData.at['leu','Percentage']}, {sampleObj.aaData.at['lys','Percentage']}, {sampleObj.aaData.at['met','Percentage']}, {sampleObj.aaData.at['phe','Percentage']}, {sampleObj.aaData.at['pro','Percentage']}, {sampleObj.aaData.at['ser','Percentage']}, {sampleObj.aaData.at['thr','Percentage']}, {sampleObj.aaData.at['trp','Percentage']}, {sampleObj.aaData.at['tyr','Percentage']}, {sampleObj.aaData.at['val','Percentage']}\n")
            file.write(f'Total AA for {domName} is {totalAA}')


def exportProteinAACounts():
    """
    Creates and writes to a new file called ProteinAACounts.csv
    File will contain Protein.proteinAACounts
    """
   #print('Beginning exportProteinAACounts')
    proteinsArray = Protein.proteinAACountsArray
    full_path = os.path.join(exportFolderName, 'ProteinAACounts.csv')
    with open(full_path, 'w') as file:
        file.write('isoform, length, totalAA, ala, arg, asn, asp, cys, gln, glu, gly, his, ile, leu, lys, met, phe, pro, ser, thr, trp, tyr, val, alaPercent, argPercent, asnPercent, aspPercent, cysPercent, gluPercent, glnPercent, glyPercent, hisPercent, ilePercent, leuPercent, lysPercent, metPercent, phePercent, proPercent, serPercent, thrPercent, trpPercent, tyrPercent, valPercent\n')
        for row in proteinsArray:
            file.write(str(row)+'\n')
        file.write('\nTotal proteins processed:'+str(len(proteinsArray)))
        file.write('\nTotal amino acids processed:'+str(Protein.AllProteinsAA))
        


def exportSampleAAComparisons():
    """
    Creates a new file with a bar chart comparing feed efficienct cows' aa profile to feed inefficient cows' profile.
    """
   #print('Beginning exportSampleAAComparisons()')
    aaSet = ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']

    inefficientData = {'ala':[],'arg':[], 'asn':[],'asp':[],'cys':[],'gln':[], 'glu':[], 'gly':[], 'his':[],'ile':[],'leu':[],'lys':[],'met':[],'phe':[],'pro':[], 'ser':[], 'thr':[], 'trp':[], 'tyr':[],'val':[]}
    efficientData = {'ala':[],'arg':[], 'asn':[],'asp':[],'cys':[],'gln':[], 'glu':[], 'gly':[], 'his':[],'ile':[],'leu':[],'lys':[],'met':[],'phe':[],'pro':[], 'ser':[], 'thr':[], 'trp':[], 'tyr':[],'val':[]}
    inefficientDataList = []
    efficientDataList = []
    df = pd.DataFrame(index = ['Rumen_1', 'Rumen_2', 'Rumen_3', 'Rumen_4', 'Rumen_5', 'Rumen_6', 'Rumen_7', 'Rumen_8', 'Rumen_9', 'Rumen_10', 'Rumen_11', 'Rumen_12', 'Rumen_13', 'Rumen_14', 'Rumen_15', 'Rumen_16', 'Rumen_17', 'Rumen_18', 'Rumen_19', 'Rumen_20', 'Rumen_21', 'Rumen_22', 'Rumen_23', 'Rumen_24'], columns = ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val', 'alaPercent','argPercent', 'asnPercent','aspPercent','cysPercent','glnPercent', 'gluPercent', 'glyPercent', 'hisPercent','ilePercent','leuPercent','lysPercent','metPercent','phePercent','proPercent', 'serPercent', 'thrPercent', 'trpPercent', 'tyrPercent','valPercent', 'Feed_Efficient?'])
    for sampObj in Sample.sampleDict.values():
        aaData = sampObj.aaData
        sampNum = int(sampObj.name[sampObj.name.rfind('-')+1:sampObj.name.find('_')])
        #  print('sample number is '+repr(sampNum))
        for aa in aaSet:
            df.at[f'Rumen_{sampNum}',aa] = aaData.at[aa,'Total']
            df.at[f'Rumen_{sampNum}',f'{aa}Percent'] = aaData.at[aa,'Percentage']
            if sampObj.feedEfficient:
                efficientData[aa].append(aaData.at[aa, 'Percentage'])
            else:
                inefficientData[aa].append(aaData.at[aa, 'Percentage'])
        df.at[f'Rumen_{sampNum}', 'Feed_Efficient?'] = sampObj.feedEfficient 
        # print('finished sample '+repr(sampNum))
        # print('df is now \n'+str(df))  
    for value in inefficientData.values():
        inefficientDataList.append(sum(value)/len(value))
    for value in efficientData.values():
        efficientDataList.append(sum(value)/len(value))
    x = np.array(range(20)) * 2
    width = 0.8
    
    fig, ax = plt.subplots(figsize = (20,6))
    rects1 = ax.bar(x - width/2, inefficientDataList, width, label = 'Feed Inefficient Samples', color = 'cornflowerblue', edgecolor = 'white')
    rects2 = ax.bar(x + width/2, efficientDataList, width, label = 'Feed Efficient Samples', color = 'green', edgecolor = 'white')
    ax.bar_label(rects1, padding = 1, fmt = '%.2f')
    ax.bar_label(rects2, padding = 1, fmt = '%.2f')

    ax.set_ylabel('Average Percentage of Microbial Protein', fontsize = 12)
    ax.set_title('Amino Acid Profiles of Feed Efficient vs Inefficient Cows', fontsize = 16)
    ax.set_xticks(x)
    for n in range(20):
        aaSet[n] = aaSet[n][0].upper()+aaSet[n][1:]
    ax.set_xticklabels(aaSet)
    ax.legend(fontsize = 13)

    plt.tight_layout()
    full_path = os.path.join(exportFolderName, 'SampleAAComparison.svg')
    plt.savefig(full_path)
    full_path = os.path.join(exportFolderName, 'SampleAAComparison.png')
    plt.savefig(full_path)
    plt.close()

    full_path = os.path.join(exportFolderName, 'SampleAAComparisonData.xlsx')
    with open(full_path, 'w') as file:
        file.write(str(df))


def exportSpeciesData():
    """
    Creates a new .csv file containing Species names as rows and amino acids as columns.
    """
   #print('Beginning exportSpeciesData()')
    full_path = os.path.join(exportFolderName, 'Species_AA_Data.csv')
    with open(full_path, 'w') as file:
        for domain, obj in Domain.domDict.items():
            file.write('\n\n--------------------------------------------------------------------\n'+domain+'\n--------------------------------------------------------------------\n\n')
            file.write('Species, Total AA, Total ala, Total arg, Total asn, Total asp, Total cys, Total gln, Total glu, Total gly, Total his, Total ile, Total leu, Total lys, Total met, Total phe, Total pro, Total ser, Total thr, Total trp, Total tyr, Total val, Percent ala, Percent arg, Percent asn, Percent asp, Percent cys, Percent gln, Percent glu, Percent gly, Percent his, Percent ile, Percent leu, Percent lys, Percent met, Percent phe, Percent pro, Percent ser, Percent thr, Percent trp, Percent tyr, Percent val\n')
            for species in obj.species.values():
                file.write(f'{species.name}, {species.totalAA}, {species.ala}, {species.arg}, {species.asn}, {species.asp}, {species.cys}, {species.gln}, {species.glu}, {species.gly}, {species.his}, {species.ile}, {species.leu}, {species.lys}, {species.met}, {species.phe}, {species.pro}, {species.ser}, {species.thr}, {species.trp}, {species.tyr}, {species.val}, {round(species.alaPercent, 2)}, {round(species.argPercent, 2)}, {round(species.asnPercent, 2)}, {round(species.aspPercent, 2)}, {round(species.cysPercent, 2)}, {round(species.glnPercent, 2)}, {round(species.gluPercent, 2)}, {round(species.glyPercent, 2)}, {round(species.hisPercent, 2)}, {round(species.ilePercent, 2)}, {round(species.leuPercent, 2)}, {round(species.lysPercent, 2)}, {round(species.metPercent, 2)}, {round(species.phePercent, 2)}, {round(species.proPercent, 2)}, {round(species.serPercent, 2)}, {round(species.thrPercent, 2)}, {round(species.trpPercent, 2)}, {round(species.tyrPercent, 2)}, {round(species.valPercent, 2)}\n')


def exportIsoformData():
    """
    Creates a new .csv file for each domain containing all isoforms within that domain, with total aa counts and percentage aa counts
    """
   #print('Beginning exportIsoformData()')
    for domName,domObj in Domain.domDict.items():
        #print(f'domName is {domName}, domObj is {domObj}')
        full_path = os.path.join(exportFolderName, f'{domName}_Isoform_Data.csv')
        with open(full_path, 'w') as file:
            #print(f'Opened {domName}_Isoform_Data.csv')
            file.write('TrinityID, Total AA, Total ala, Total arg, Total asn, Total asp, Total cys, Total gln, Total glu, Total gly, Total his, Total ile, Total leu, Total lys, Total met, Total phe, Total pro, Total ser, Total thr, Total trp, Total tyr, Total val, Percent ala, Percent arg, Percent asn, Percent asp, Percent cys, Percent gln, Percent glu, Percent gly, Percent his, Percent ile, Percent leu, Percent lys, Percent met, Percent phe, Percent pro, Percent ser, Percent thr, Percent trp, Percent tyr, Percent val\n')
            for protObj in domObj.isoforms.values():
                file.write(f'{protObj.id}, {protObj.totalAA}, {protObj.ala}, {protObj.arg}, {protObj.asn}, {protObj.asp}, {protObj.cys}, {protObj.gln}, {protObj.glu}, {protObj.gly}, {protObj.his}, {protObj.ile}, {protObj.leu}, {protObj.lys}, {protObj.met}, {protObj.phe}, {protObj.pro}, {protObj.ser}, {protObj.thr}, {protObj.trp}, {protObj.tyr}, {protObj.val}, {round(protObj.alaPercent, 2)}, {round(protObj.argPercent, 2)}, {round(protObj.asnPercent, 2)}, {round(protObj.aspPercent, 2)}, {round(protObj.cysPercent, 2)}, {round(protObj.glnPercent, 2)}, {round(protObj.gluPercent, 2)}, {round(protObj.glyPercent, 2)}, {round(protObj.hisPercent, 2)}, {round(protObj.ilePercent, 2)}, {round(protObj.leuPercent, 2)}, {round(protObj.lysPercent, 2)}, {round(protObj.metPercent, 2)}, {round(protObj.phePercent, 2)}, {round(protObj.proPercent, 2)}, {round(protObj.serPercent, 2)}, {round(protObj.thrPercent, 2)}, {round(protObj.trpPercent, 2)}, {round(protObj.tyrPercent, 2)}, {round(protObj.valPercent, 2)}\n')
                #print(f'Wrote {protObj.id} to file')


def exportSamplePieCharts():
    """
    Creates a new file named 'Sample_Pie_Charts.png' that contains a pie chart for each sample.
    
    Creates a new file named 'Averaged_Pie_Chart.png' that contains a single pie chart of the average values of the samples.
    
    Each pie chart represents the total number of amino acids by domain.
    """
   #print('Beginning exportSamplePieCharts()')
    i = 1
    plt.figure(figsize=(14,len(Sample.sampleDict)*1.8))
    plt.subplots_adjust(hspace = 1.2)
    avg_values = [0,0,0,0,0]
    labels = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'Unknown']
    for sampleObj in Sample.sampleDict.values():
        archaeaSum = sampleObj.aaData['Total_Archaea'].sum()
        bacteriaSum = sampleObj.aaData['Total_Bacteria'].sum()
        eukaryotaSum = sampleObj.aaData['Total_Eukaryota'].sum()
        virusesSum = sampleObj.aaData['Total_Viruses'].sum()
        unknownSum = sampleObj.aaData['Total_Unknown'].sum()
        values = np.array([archaeaSum,bacteriaSum,eukaryotaSum,virusesSum,unknownSum])
        avg_values[0] += archaeaSum
        avg_values[1] += bacteriaSum
        avg_values[2] += eukaryotaSum
        avg_values[3] += virusesSum
        avg_values[4] += unknownSum

        #print('values are:\n'+str(values))
        #print('avg_values are:\n'+str(avg_values))
        if values.sum() != 0:
            plt.subplot(-(-len(Sample.sampleDict)//2), 2, i)
            plt.pie(values, radius = 1.5, autopct = '%1.1f%%')
            plt.title(sampleObj.name, pad = 25)
            plt.legend(labels, bbox_to_anchor = (1.2,0.7))
        i += 1
    plt.tight_layout()
    full_path = os.path.join(exportFolderName, 'Sample_Pie_Charts.png')
    plt.savefig(full_path)
    full_path = os.path.join(exportFolderName, 'Sample_Pie_Charts.svg')
    plt.savefig(full_path)
    plt.close()

    plt.figure(figsize = (10,10))
    plt.pie(np.array(avg_values), radius = 1, autopct = '%1.1f%%')
    plt.title('Average Sample Amino Acid Sources')
    plt.legend(labels)
    plt.tight_layout()
    full_path = os.path.join(exportFolderName, 'Average_Sample_Pie_Chart.svg')
    plt.savefig(full_path)
    full_path = os.path.join(exportFolderName, 'Average_Sample_Pie_Chart.png')
    plt.savefig(full_path)



def exportAAHistograms(target_aa):
    """
    Creates a new file for each sample named '{Sample}_AAs_by_Isoforms.png' where {Sample} is the name of a sample.
    Each file contains 20 histograms, one for each amino acid, plotting number of isoforms on the y axis by percent abundance of amino acid
    on the x axis.
    For target_aa, creates a new file called {aa}Histogram, with an independent image of the histogram for that aa% in proteins across all domains.

    target_aa is a string of one of the 20 three-letter abbreviations for amino acids, all lowercase.
    """
   #print('Beginning exportSampleHistograms()')
    aaSet = ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']
    
    aaHistData = []
    for domObj in Domain.domDict.values():
        fig, ax = plt.subplots(figsize = (12,20))
        i = 1
        data = {'ala':[],'arg':[], 'asn':[],'asp':[],'cys':[], 'gln':[], 'glu':[], 'gly':[], 'his':[],'ile':[],'leu':[],'lys':[],'met':[],'phe':[],'pro':[], 'ser':[], 'thr':[], 'trp':[], 'tyr':[],'val':[]}
        for aa in aaSet:
            plt.subplot(10, 2, i)
            for protObj in domObj.isoforms.values():
                data[aa].append(getattr(protObj, f'{aa}Percent'))
                if aa == target_aa:
                    aaHistData.append(getattr(protObj, f'{aa}Percent'))
            narrowData = np.array(data[aa])
            counts, bins, bars = plt.hist(narrowData, bins = 20, edgecolor = 'white')
            plt.bar_label(bars, padding =5)
            plt.title(aa[0].upper()+aa[1:])
            plt.ylabel('# of Isoforms')
            plt.xlabel(f'% {aa} in Isoform')
            
            i += 1
        fig.suptitle(domObj.name)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 1.5, wspace = 0.2)

        full_path = os.path.join(exportFolderName, f'{domObj.name}_AA_Histograms.png')
        plt.savefig(full_path)
        plt.close()
    plt.figure(figsize = (20,8))
    counts, bins, bars = plt.hist(aaHistData, bins = 30, edgecolor = 'white')
    plt.bar_label(bars, padding = 10)
    plt.title(target_aa.upper(), fontsize = 24)
    plt.ylabel('# of Isoforms', fontsize = 18)
    plt.xlabel(f'Percent {target_aa} in Isoform', fontsize = 18)
    plt.tight_layout()

    full_path = os.path.join(exportFolderName, f'{target_aa.upper()}IndependentHistogram.svg')
    plt.savefig(full_path)
    full_path = os.path.join(exportFolderName, f'{target_aa.upper()}IndependentHistogram.png')
    plt.savefig(full_path)
    plt.close()

def exportSpeciesAAHistogram(target_aa):
    """
    Creates a new file called {AA}SpeciesHistogram.png with a single histogram plotting number of species on the y axis
    and percent of target_aa on the x axis.

    target_aa is a lowercase string of one of the 20 three-letter amino acid abbreviations.
    """
   #print('Beginning exportSpeciesAAHistogram')
    assert isinstance(target_aa, str)

    data = []
    full_path = os.path.join(exportFolderName, f'{target_aa.upper()}SpeciesHistogramDATA.csv')
    with open(full_path, 'w') as file:
        file.write('Species,Percentage\n')
        for specObj in Species.speciesDict.values():
            if specObj.totalAA > 0:
                percent = getattr(specObj, target_aa)/specObj.totalAA * 100
                if percent > 0:
                    data.append(percent)
                    file.write(specObj.name+','+str(percent)+'\n')
        
    plt.figure(figsize = (20,8))
    counts, bins, bars = plt.hist(data, bins = 30, edgecolor = 'white')
    plt.bar_label(bars, padding = 10)
    plt.title(target_aa.upper(), fontsize = 24)
    plt.ylabel('Microbial Species in the rumen, n', fontsize = 18)
    plt.xlabel(f'{target_aa} in their Protein Composition, %', fontsize = 18)
    plt.tight_layout()

    full_path = os.path.join(exportFolderName, f'{target_aa.upper()}SpeciesHistogram.svg')
    plt.savefig(full_path)
    full_path = os.path.join(exportFolderName, f'{target_aa.upper()}SpeciesHistogram.png')
    plt.savefig(full_path)
    plt.close()
    

def exportSelectAAData(aa):
    """
    Creates two new files for each domain, titled {Domain}Lowest{aa}Proteins.txt and {Domain}Highest{aa}Proteins.txt.
    The Lowest file contains the lowest 10% proteins in the domain for aa and the Highest file contains the top 10%.

    aa must be a string of an amino acid abbreviation. Any of ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']
    """
   #print('Beginning exportSelectAAData()')
    assert isinstance(aa, str)
    assert len(aa) > 0
    assert aa in ['ala','arg', 'asn','asp','cys','gln', 'glu', 'gly', 'his','ile','leu','lys','met','phe','pro', 'ser', 'thr', 'trp', 'tyr','val']

    aaPercents = []
    for protObj in Protein.proteinDict.values():
        aaPercents.append(getattr(protObj, f'{aa}Percent'))
    aaPercents = np.array(aaPercents)
    cutoff1 = np.percentile(aaPercents, 1)
    cutoff99 = np.percentile(aaPercents, 99)
    bottom10 = {}
    upper10 = {}
    for protObj in Protein.proteinDict.values():
        protPercent = round(getattr(protObj, f'{aa}Percent'),2)
        if protPercent <= cutoff1:
            bottom10[protObj.id] = str(protPercent)
        elif protPercent >= cutoff99:
            upper10[protObj.id] = str(protPercent)

    full_path = os.path.join(exportFolderName, f'Lowest{aa.upper()}Proteins.txt')
    with open(full_path, 'w') as file:
        file.write('Isoform,Percentage\n')
        for id, percent in bottom10.items():
            file.write(id+','+percent+'\n')
    full_path = os.path.join(exportFolderName, f'Highest{aa.upper()}Proteins.txt')
    with open(full_path, 'w') as file:
        file.write('Isoform,Percentage\n')
        for id, percent in upper10.items():
            file.write(id+','+percent+'\n')






