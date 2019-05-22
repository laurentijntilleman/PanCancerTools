#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 14:15:02 2019

@author: Laurentijn Tilleman
"""

path = 'data/'

import pandas as pd
import requests
import cosmic_lookup_table

# load all mutations from Ellrott et al.
mutations = pd.read_csv(path+'Mutation.CTAT.3D.Scores.txt',sep='\t')

# select driver mutations
driverMutations = mutations[mutations['New_Linear (cancer-focused) flag']+
                            mutations['New_Linear (functional) flag']+
                            mutations['3D mutational hotspot flag']>1]

# cosmic lookup table
LOOKUP_TABLE = cosmic_lookup_table.CosmicLookup(path+"cosmic_lookup_table.tsv")
# Ensembl server
server = "http://grch37.rest.ensembl.org"

driverMutations['chr'] = ''
driverMutations['start'] = 0
driverMutations['end'] = 0
driverMutations['ref'] = ''
driverMutations['alt'] = ''
driverMutations['build'] = ''

for index in driverMutations.index:
    gene = driverMutations.loc[index,'gene']
    alteration = driverMutations.loc[index,'protein_change']
    matches = LOOKUP_TABLE.get_entries(gene,alteration)
    if len(matches) > 0:
        driverMutations.loc[index,'chr'] = matches[0]['chrom']
        driverMutations.loc[index,'start'] = matches[0]['start']
        driverMutations.loc[index,'end'] = matches[0]['end']
        driverMutations.loc[index,'ref'] = matches[0]['ref']
        driverMutations.loc[index,'alt'] = matches[0]['alt']
        driverMutations.loc[index,'build'] = str(int(matches[0]['build']))
    else:
        ext = "/vep/human/hgvs/{}:{}?".format(gene,alteration)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if r.ok:
            decoded = r.json()
            driverMutations.loc[index,'chr'] = decoded[0]['seq_region_name']
            driverMutations.loc[index,'start'] = decoded[0]['start']
            driverMutations.loc[index,'end'] = decoded[0]['end']
            driverMutations.loc[index,'ref'] = decoded[0]['allele_string'].split('/')[0]
            driverMutations.loc[index,'alt'] = decoded[0]['allele_string'].split('/')[1]
            driverMutations.loc[index,'build'] = decoded[0]['assembly_name']

driverMutations.to_csv(path+'driverMutation.csv')
