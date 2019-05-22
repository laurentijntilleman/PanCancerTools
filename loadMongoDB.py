#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 15:30:04 2019

@author: Laurentijn Tilleman
"""

import json
from pymongo import MongoClient

path = 'g2p-aggregator/harvester/'
datafiles = ['cgi.json','civic.json','oncokb.json','depo.json','jax.json']

client = MongoClient()
db = client.cancerPGx

db.PGxAll.remove()

for data in datafiles:
    with open(path+data) as fh:
        line = fh.readline()
        while line:
            db.PGxAll.insert(json.loads(line))
            line = fh.readline()
