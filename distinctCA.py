#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:25:19 2019

@author: Laurentijn Tilleman
"""

from pymongo import MongoClient
import pandas as pd
import re
import CrossMap

# connect to mongoDB
client = MongoClient()
db = client.cancerPGx

# responce types
sensitive = ['Responsive','Sensitivity/Response','sensitive']
resistance = ['Reduced Sensitivity','resistant','No Responsive','no benefit']
toxic = ["Increased Toxicity","Increased Toxicity (Haemolytic Anemia)",
         "Increased Toxicity (Myelosupression)",
         "Increased Toxicity (Ototoxicity)",
         "Increased Toxicity (Hyperbilirubinemia)",
         "Adverse Response"]
other = ['Positive','Likely Pathogenic','Better Outcome','Negative',None,
         "Pathogenic",'Uncertain Significance','N/A','Loss of Function',
         'Gain of Function',"Neomorphic"]

# initialization
CA = dict()
variantNotComplet = 0

# used responce types
responceTypes = ['sensitive','resistance','toxic']

# list of drugs correct split
comma = list(set(pd.read_csv('data/comma2.txt',sep='\t',header=None)[0]))

# filtering duplicated entries
for responceType in responceTypes:
    CA[responceType] = dict()
    results = db.PGxAll.aggregate([{'$unwind':'$genes'},
            {"$match":{"association.drug_labels":{"$exists":True}}},
            {"$match":{"association.response_type":{"$in":eval(responceType)}}}
            ])
    for result in results:
        if result['features'] != [] and \
            all (k in resultK for resultK in result['features'] for k in ('chromosome','start','end')) \
            and all ( resultK[k] != None  for resultK in result['features'] for k in ('chromosome','start','end')):
            if result['genes'] not in CA[responceType]:
                CA[responceType][result['genes']] = dict()
            if ',' in result['association']['drug_labels'] and any(x.isdigit() for x in result['association']['drug_labels']):
                drugs=[commaString for commaString in comma if commaString in result['association']['drug_labels']]
            else:
                drugs = result['association']['drug_labels'].split(',')
                for drug in drugs:
                    drug = drug.lower()
                    if drug not in CA[responceType][result['genes']]:
                        CA[responceType][result['genes']][drug] = dict()
                    featureDescription = '_'.join(sorted(set([re.sub(' +',' ',feature['description']) for feature in result['features']])))
                    if featureDescription not in CA[responceType][result['genes']][drug]:
                        CA[responceType][result['genes']][drug][featureDescription] = {'evidence':[result['association']['evidence_label']],
                          'features':result['features'],
                          'phenotypes':set([x['description'] for x in result['association']['phenotypes']]),
                          'source':[result['source']]}
                    else:
                        if result['source'] not in CA[responceType][result['genes']][drug][featureDescription]['source']:
                            CA[responceType][result['genes']][drug][featureDescription]['source'].append(result['source'])
                        CA[responceType][result['genes']][drug][featureDescription]['evidence'].append(result['association']['evidence_label'])
                        phenotypes = CA[responceType][result['genes']][drug][featureDescription]['phenotypes'] | set([x['description'] for x in result['association']['phenotypes']])
                        CA[responceType][result['genes']][drug][featureDescription]['phenotypes']=phenotypes
        else:
            variantNotComplet +=1

# load database
data = []
for responceType in responceTypes:
    for gene in CA[responceType]:
        for drug in CA[responceType][gene]:
            for feature in CA[responceType][gene][drug]:
               # if CA[responceType][gene][drug][feature]['conflicting'] == 0:
                    featureDict = CA[responceType][gene][drug][feature]
                    evidence = sorted(featureDict['evidence'])[0]
                    data.append({'gene':gene,'drug':drug,
                        'responceType':responceType,
                        'features':featureDict['features'],
                        'evidenceLabel':evidence,
                        'source':featureDict['source'],
                        'phenotypes':[item for sublist in featureDict['phenotypes'] for item in sublist.split(';')]})
db.PGxDD.remove()
db.PGxDD.insert_many(data)

# chain file to convert the genomic cordinates with CrossMap
chain_file = 'data/hg38ToHg19.over.chain'
(mapTree,targetChromSizes, sourceChromSizes) = CrossMap.read_chain_file(chain_file)

# UCSC exon coordinates
exons = pd.read_csv('data/UCSC_exons_modif_canonical.bed',header=None,sep='\t')
intron = pd.read_csv('data/UCSC_introns_modif_canonical.bed',header=None,sep='\t')
exons3 = pd.read_csv('data/UCSC_3_exons_modif_canonical.bed',header=None,sep='\t')

# convert the genomic coordinates from the AVENIO panels
avenioTargeted = pd.read_excel('data/Targeted_Panel_Regions_Information_Tumor_08374317001.xlsx',header=None)
for line in avenioTargeted.index:
    output = CrossMap.map_coordinates(mapTree, avenioTargeted.loc[line][0], avenioTargeted.loc[line][1], avenioTargeted.loc[line][2])
    avenioTargeted.loc[line,0] = output[-1][0]
    avenioTargeted.loc[line,1] = output[-1][1]
    avenioTargeted.loc[line,2] = output[-1][2]

avenioExpanded = pd.read_excel('data/Expanded_Panel_Regions_Information_Tumor_08374325001.xlsx',header=None)
for line in avenioExpanded.index:
    output = CrossMap.map_coordinates(mapTree, avenioExpanded.loc[line][0], avenioExpanded.loc[line][1], avenioExpanded.loc[line][2])
    avenioExpanded.loc[line,0] = output[-1][0]
    avenioExpanded.loc[line,1] = output[-1][1]
    avenioExpanded.loc[line,2] = output[-1][2]

avenioSurveillence = pd.read_excel('data/Surveillance_Panel_Regions_Information_Tumor_08374333001.xlsx',header=None)
for line in avenioSurveillence.index:
    output = CrossMap.map_coordinates(mapTree, avenioSurveillence.loc[line][0], avenioSurveillence.loc[line][1], avenioSurveillence.loc[line][2])
    avenioSurveillence.loc[line,0] = output[-1][0]
    avenioSurveillence.loc[line,1] = output[-1][1]
    avenioSurveillence.loc[line,2] = output[-1][2]

# load FoundationMedicine genes panels
oneL = set(pd.read_csv('data/FoundationOne.txt', sep=' *, *',header=None)[0])
oneCDx = set()
with open('data/FoundationOneCDx.txt','r') as fh:
    for line in fh.readlines():
        for gene in re.sub("\(.*?\)",'',line).replace('*','').split():
            oneCDx.add(gene)
oneCDxIntron = pd.read_csv('data/FoundationOneCDxIntron.txt',sep='\t',header=None)
oneCDx3Exon = pd.read_csv('data/FoundationOneCDx3Exon.txt',sep='\t',header=None)
oneL_Intron = pd.read_csv('data/FoundationOneIntron.txt',sep='\t',header=None)

# lookup genomic coordinates of the selected genes
foundationOneCDx = pd.DataFrame()
for line in oneCDx:
    foundationOneCDx = foundationOneCDx.append(exons[exons[3] == line.strip()])
for line in oneCDxIntron.index:
    foundationOneCDx = foundationOneCDx.append(intron[intron[3]==oneCDxIntron.loc[line][0]][intron[4]==int(oneCDxIntron.loc[line][1])-1])
for line in oneCDx3Exon.index:
    foundationOneCDx = foundationOneCDx.append(exons3[exons3[3] == oneCDx3Exon.loc[line][0]])
foundationOneLiquid = pd.DataFrame()
for line in oneL:
    foundationOneLiquid = foundationOneLiquid.append(exons[exons[3] == line.strip()])
for line in oneL_Intron.index:
    foundationOneLiquid = foundationOneLiquid.append(intron[intron[3]==oneL_Intron.loc[line][0]][intron[4]==int(oneL_Intron.loc[line][1])-1])

# load all pannels
panels = [pd.read_csv('data/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.bed',header=None,sep='\t'),
          pd.read_csv('data/PanCancerV2.4.target(hg19).bed',header=None,sep='\t'),
          pd.read_csv('data/TST500C_manifest.bed',header=None,sep='\t'),
          pd.read_csv('data/CCP.20131001.submitted.bed',header=None,sep='\t',skiprows=1),
          pd.read_csv('data/DHS-3501Z.covered-250bp.bed',header=None,sep='\t',skiprows=1),
          avenioTargeted,
          avenioExpanded,
          avenioSurveillence,
          foundationOneCDx[~foundationOneCDx.index.duplicated(keep='first')],
          foundationOneLiquid
          ]

# pan cancer datasets
panCancerMutations = pd.read_csv('data/driverMutation.csv')
panCancerFusions = pd.read_excel('data/1-s2.0-S2211124718303954-mmc2.xlsx',sheetname=1,skiprows=1)

# define the index
dataSets = ['sensitive A','sensitive B','sensitive C','sensitive D',
            'resistance A','resistance B','resistance C','resistance D',
            'toxic A','toxic B','toxic C','toxic D',
            'PanCancer Mutations','PanCancer Fusion']

# define the columns
columns = ['Total','Exome IDT','IDT2','Illumina','AmpliSeq','Qiagen',
           'AVENIO Targeted','AVENIO Expanded','AVENIO Surveillance',
           'FoundationOne CDx','FoundationOne Liquid']

# initialize empty dataframe
dataTable = pd.DataFrame(0,index = dataSets +
                                    list(set([x for y in set(panCancerMutations['CODE']) for x in y.split(',')])),
                         columns = columns)

# query all pharmacogenomic interactions
data = [x for x in db.PGxDD.find({})]

# count pharmacogenomic interactions
for d in data:
    columnN = 1
    dataTable['Total'][d['responceType']+' '+d['evidenceLabel'][0]] +=1
    plusCount = False
    for panel in panels:
        count = 0
        column = dataTable.columns[columnN]
        if len(d['features']) == 1:
            if 'end' in d['features'][0] and 'start' in d['features'][0]:
                feature = d['features'][0]
                try:
                    feature['chromosome'] = int(float(feature['chromosome']))
                except ValueError:
                    pass
                if d['features'][0]['end'] == d['features'][0]['start']:
                    if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) != 0:
                        count += 1
                else:
                    if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) != 0:
                        count += 1
                    elif 'provenance_rule' in feature:
                        if feature['provenance_rule'] in ['is_amplification','is_deletion','is_loss']:
                            if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                count += 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                count += 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                count += 1
                    elif 'name' in feature:
                        if feature['name'] in ['loss']:
                            if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                count += 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                count += 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                count += 1
        else:
            count2 = True
            for feature in d['features']:
                if 'end' in feature and 'start' in feature:
                    try:
                        feature['chromosome'] = int(float(feature['chromosome']))
                    except ValueError:
                        pass
                    if feature['end'] != feature['start']:
                        if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) != 0:
                            pass
                        elif 'biomarker_type' in feature:
                            if feature['biomarker_type'] in ['amp','del']:
                                if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                    pass
                                elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                    pass
                                elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                    pass
                                else:
                                    count2 = False
                            else:
                                count2 = False
                        elif 'name' in feature:
                            if feature['name'] in ['loss']:
                                if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                    pass
                                elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                    pass
                                elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                    pass
                                else:
                                    count2 = False
                        else:
                            count2 = False
                    elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) == 0:
                        count2 = False
                else:
                    count2 = False
            if count2:
                count += 1
        if count > 0:
            dataTable[column][d['responceType']+' '+d['evidenceLabel'][0]] +=1
        columnN +=1

# count Pan-Cancer fusions
countFusions = 0
for line in panCancerFusions.index:
    chromosome = panCancerFusions.loc[line,'Breakpoint1'].split(':')[0]
    start = int(panCancerFusions.loc[line,'Breakpoint1'].split(':')[1])
    end = int(panCancerFusions.loc[line,'Breakpoint1'].split(':')[1])
    strand = panCancerFusions.loc[line,'Breakpoint1'].split(':')[2]
    output = CrossMap.map_coordinates(mapTree, chromosome, start, end, strand)
    if output == None:
        countFusions += 1
        output1 = (chromosome,start,end,strand)
    else:
        output1 = output[-1]
    chromosome = panCancerFusions.loc[line,'Breakpoint2'].split(':')[0]
    start = int(panCancerFusions.loc[line,'Breakpoint2'].split(':')[1])
    end = int(panCancerFusions.loc[line,'Breakpoint2'].split(':')[1])
    strand = panCancerFusions.loc[line,'Breakpoint2'].split(':')[2]
    output = CrossMap.map_coordinates(mapTree, chromosome, start, end, strand)
    if output == None:
        countFusions += 1
        output2 = (chromosome,start,end,strand)
    else:
        output2 = output[-1]
    dataTable['Total']['PanCancer Fusion'] +=1
    columnN = 1
    for panel in panels:
        column = dataTable.columns[columnN]
        if len(panel[panel[0] == output1[0]][panel[1] <= output1[1]][panel[2] >= output1[2]]) != 0:
            if len(panel[panel[0] == output2[0]][panel[1] <= output2[1]][panel[2] >= output2[2]]) != 0:
                dataTable[column]['PanCancer Fusion'] += 1
        columnN +=1

# count Pan-Cancer mutations
for line in panCancerMutations.index:
    dataTable['Total']['PanCancer Mutations'] +=1
    for cancer in panCancerMutations.loc[line,'CODE'].split(','):
        dataTable['Total'][cancer] += 1
    columnN = 1
    try:
        panCancerMutations.loc[line,'chr'] = int(float(panCancerMutations.loc[line,'chr']))
    except ValueError:
        pass
    for panel in panels:
        column = dataTable.columns[columnN]
        if len(panel[panel[0] == 'chr{}'.format(panCancerMutations.loc[line,'chr'])][panel[1] <= panCancerMutations.loc[line,'start']][panel[2] >= panCancerMutations.loc[line,'end']]) != 0:
            dataTable[column]['PanCancer Mutations'] += 1
            for cancer in panCancerMutations.loc[line,'CODE'].split(','):
                dataTable[column][cancer] += 1
        columnN +=1

dataTable.to_csv('results/dataTable.csv')

## Variants
# load all pharmacogenomic interactions
data = db.PGxDD.find({})

count = 0

source = dict()
variants = dict()

variantsDict = dict()

for d in data:
    try:
        source['_'.join(sorted(d['source']))] +=1
    except KeyError:
        source['_'.join(sorted(d['source']))] =1
    for feature in d['features']:
        variantsDict[re.sub(' +',' ',feature['description'])] = feature
        try:
            for s in d['source']:
                variants[re.sub(' +',' ',feature['description'])].add(s)
        except KeyError:
                variants[re.sub(' +',' ',feature['description'])] = set()
                for s in d['source']:
                    variants[re.sub(' +',' ',feature['description'])].add(s)

# insert to database
db.variants.remove()
db.variants.insert_many([variantsDict[x] for x in variantsDict])

# variants per panel
panelVariants = pd.DataFrame(data = 0, columns = columns, index= ['variants'])
UpSet_panelVariants = pd.DataFrame(columns = columns[1:])
variantCount = 0
variantCountN = 0

for variant in variantsDict:
    if all (k in variantsDict[variant] for k in ('chromosome','start','end')):
        feature = variantsDict[variant]
        if all( feature[k] != 'None' for k in ('chromosome','start','end')):
            try:
                feature['chromosome'] = int(float(feature['chromosome']))
            except ValueError:
                pass
            panelVariants.loc['variants','Total'] +=1
            panelCol = 1
            upSetList = [0] * 11
            variantCount += 1
            for panel in panels:
                if feature['end'] == feature['start']:
                    if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) != 0:
                        panelVariants.loc['variants',columns[panelCol]] +=1
                        upSetList[panelCol-1] = 1
                else:
                    if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] <= feature['start']][panel[2] >= feature['end']]) != 0:
                        panelVariants.loc['variants',columns[panelCol]] +=1
                        upSetList[panelCol-1] = 1
                    elif 'provenance_rule' in feature:
                        if feature['provenance_rule'] in ['is_amplification','is_deletion','is_loss']:
                            if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                    elif 'name' in feature:
                        if feature['name'] in ['loss']:
                            if len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[1] >= feature['start']][panel[1] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                            elif len(panel[panel[0] == 'chr{}'.format(feature['chromosome'])][panel[2] >= feature['start']][panel[2] <= feature['end']]) != 0:
                                panelVariants.loc['variants',columns[panelCol]] +=1
                                upSetList[panelCol-1] = 1
                panelCol +=1
            exec('UpSet_{}.loc["{}"] = {}'.format('panelVariants',variantCount,upSetList))

UpSet_panelVariants.to_csv("results/UpSet_panelVariants.csv")       
panelVariants.to_csv('results/panelVariants.csv')

# count variants in database
numVar = dict()
for variant in variants:
    try:
        numVar[len(variants[variant])] += 1
    except KeyError:
        numVar[len(variants[variant])] = 1

## sources of knowledge bases

sourceTable = pd.DataFrame(columns = ['cgi','civic','depo','jax','oncokb'])

for x in db.PGxDD.find({}):
    sourceTable.loc[x['_id']] = [1 if y in x['source'] else 0 for y in sourceTable.columns]

sourceTable.to_csv('results/sourceTable.csv')

## Overlap meta-knowledge base and Pan-Cancer driver mutations

overlap = 0

def reforme(x):
    try:
        return int(float(x))
    except ValueError:
        return x

panCancerMutations['chr'] = [reforme(x) for x in panCancerMutations['chr']]
panCancerMutations['start'] = [int(x) for x in panCancerMutations['start']]
panCancerMutations['end'] = [int(x) for x in panCancerMutations['end']]

for variant in variantsDict:
    if all (k in variantsDict[variant] for k in ('chromosome','start','end')):
        feature = variantsDict[variant]
        if all( feature[k] != 'None' for k in ('chromosome','start','end')) \
            and all ( feature[k] != None for k in ('chromosome','start','end')):
            try:
                feature['chromosome'] = int(float(feature['chromosome']))
            except ValueError:
                pass
            if len(panCancerMutations[panCancerMutations['chr'] == \
                                      feature['chromosome']][panCancerMutations['start'] \
                                             <= int(feature['start'])][panCancerMutations['end'] \
                                            >= int(feature['end'])]) != 0:
                overlap += 1
