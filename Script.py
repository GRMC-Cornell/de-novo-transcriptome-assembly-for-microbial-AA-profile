import Funcs
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 10000)
pd.set_option('display.max_colwidth', None)
Funcs.createDomains()



taxonomyFile =  input("Enter the Taxonomy File name. Make sure it is in the same folder as this script. \n") 
proteinFile = input("Enter the Protein Sequence File name. Make sure it is in the same folder as this script.\n") 
abundanceFile = input('Enter the abundance File name. Make sure it is in the same folder as this script. \n')
target_aa = input('Choose which amino acid to generate graphs. Enter its three-letter abbreviation: \n')
target_aa = target_aa.lower()

print('Beginning taxonomy file extraction and analysis')
# Read, parse taxonomy file. Create Species and Domain objects.
taxonomyTable = Funcs.extractBLAST(taxonomyFile)
catalogue = Funcs.taxParse(taxonomyTable)
Funcs.assignSpecies(catalogue)

print('Beginning protein sequence file extraction and analysis')
# Read, parse dataset. Create Protein objects, assign them to Species and Domain objects.
rawDataset = Funcs.extractDataset(proteinFile)
protMatch = Funcs.transformDataset(rawDataset)
Funcs.assignIsoform(protMatch, catalogue)

print('Beginning abundance file extraction and analysis')
# Read, parse abundance file. Create Sample objects, apply abundance values to proteins for each sample.
abundanceData = Funcs.extractAbundance(abundanceFile)
Funcs.createSamples(abundanceData)
Funcs.applySampleMetadata()
Funcs.differentiateAbundance(abundanceData)
Funcs.computeRawAA(protMatch)
Funcs.computePercentAA()
Funcs.computeSpeciesAACounts()

# Export results
Funcs.createExportFolder()
Funcs.exportProteinAACounts()
Funcs.exportSampleData()
Funcs.exportSpeciesData()
Funcs.exportIsoformData()
Funcs.exportSampleAAComparisons()
Funcs.exportSamplePieCharts()
Funcs.exportAAHistograms(target_aa)
Funcs.exportSpeciesAAHistogram(target_aa)
Funcs.exportSelectAAData(target_aa)
print('\n\n\nArchaea species are: \n'+ str(Funcs.Domain.domDict['Archaea'].species.keys()))
#print('Archaea isoforms are: \n'+ str(Funcs.Domain.domDict['Archaea'].isoforms))

print('\n\n\nBacteria species are: \n'+ str(Funcs.Domain.domDict['Bacteria'].species.keys()))
#print('Bacteria isoforms are: \n'+ str(Funcs.Domain.domDict['Bacteria'].isoforms))

print('\n\n\nEukaryota species are: \n'+ str(Funcs.Domain.domDict['Eukaryota'].species.keys()))
#print('Eukaryota isoforms are: \n'+ str(Funcs.Domain.domDict['Eukaryota'].isoforms))

print('\n\n\nViruses species are: \n'+ str(Funcs.Domain.domDict['Viruses'].species.keys()))
#print('Viruses isoforms are: \n'+ str(Funcs.Domain.domDict['Viruses'].isoforms))

print('\n\n\nUnknown species are: \n'+ str(Funcs.Domain.domDict['Unknown'].species.keys()))
#print('Unknown isoforms are: \n'+ str(Funcs.Domain.domDict['Unknown'].isoforms))

print('\n\n\nSamples are: \n'+str(Funcs.Sample.sampleDict.keys()))


print('\n\nJob Done. New files created.')

