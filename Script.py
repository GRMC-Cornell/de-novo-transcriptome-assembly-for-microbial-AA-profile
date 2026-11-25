import Funcs
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 10000)
pd.set_option('display.max_colwidth', None)




taxonomyFile =  input("Enter the Taxonomy File name. Make sure it is in the same folder as this script. \n") 
proteinFile = input("Enter the Protein Sequence File name. Make sure it is in the same folder as this script.\n") 
abundanceFile = input('Enter the transcript TPM matrix filename (tab-delimited; transcript IDs in the first column, samples in the remaining columns):\n').strip()
target_aa = input('Choose which amino acid to generate graphs. Enter its three-letter abbreviation: \n')
target_aa = target_aa.lower()

#### To automate the script further, de-comment the following lines ####
#### after entering your file names into the respective sections.   ####
#### Then feel free to comment out the lines above this.            ####

# taxonomyFile = 'Enter your file path here'
# proteinFile = 'Enter your file path here'
# abundanceFile = 'Enter your file path here'
# target_aa = 'Enter the 3 letter abbreviation here'

# First, create the folder for exports.
Funcs.createExportFolder()

# Read, parse taxonomy file. Create Species and Domain objects.
print('Beginning taxonomy file extraction and analysis')
Funcs.createDomains()
taxonomyTable = Funcs.extractBLAST(taxonomyFile)
catalogue = Funcs.taxParse(taxonomyTable)
Funcs.assignSpecies(catalogue)

# Read, parse dataset. Create Protein objects, assign them to Species and Domain objects.
print('Beginning protein sequence file extraction and analysis')
rawDataset = Funcs.extractDataset(proteinFile)
protMatch = Funcs.transformDataset(rawDataset)
Funcs.assignIsoform(protMatch, catalogue)

# Read, parse abundance file. Create Sample objects, apply abundance values to proteins for each sample.
print('Beginning abundance file extraction and analysis')
abundanceData = Funcs.extractAbundance(abundanceFile)
Funcs.createSamples(abundanceData)
Funcs.applySampleMetadata()
Funcs.differentiateAbundance(abundanceData)

# Apply abundance multipliers to amino acid counts
print('Combining abundance data with raw amino acid counts')
Funcs.computeRawAA(protMatch)
Funcs.computePercentAA()
Funcs.computeSpeciesAACounts()

# Export results
print('Creating Files and Exporting Results')
Funcs.exportProteinAACounts()
Funcs.exportSampleData()
Funcs.exportSpeciesData()
Funcs.exportIsoformData()
Funcs.exportSampleAAComparisons()
Funcs.exportSamplePieCharts()
Funcs.exportAAHistograms(target_aa)
Funcs.exportSpeciesAAHistogram(target_aa)
Funcs.exportSelectAAData(target_aa)


print('\n\nJob Done. New files created in GRMC_Exports.')

