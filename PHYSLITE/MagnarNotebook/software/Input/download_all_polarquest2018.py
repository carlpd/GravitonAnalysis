#!/usr/bin/env python
# Script to download all .nc files from a THREDDS catalog directory
# Written by Sage 4/5/2016, revised 5/31/2018

import wget
from datetime import timedelta, date
import requests
import sys
from bs4 import BeautifulSoup
from urllib import *

while True:
    g = input("Specify file format to download (1: csv, 2: root, 3: both)") 
    try:
        filef = int(g)
        break
    except:
        print("Please specify an integer in [1-3]")
        
        
folder_list = []
file_list = []
def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

start_date = date(2018, 5, 28)
end_date = date(2021, 8, 13)
for single_date in daterange(start_date, end_date):
    #print(single_date.strftime("%Y-%m-%d"))
    folder_list.append(single_date.strftime("%Y-%m-%d"))

for i in range(len(folder_list)):
    
    url = "https://iatw.cnaf.infn.it/eee/monitor/dqmreport2/dqmreport2/POLA-03/"+folder_list[i]

    print("\nDownloading file %i/%i (now: %s)" %(i,len(folder_list),folder_list[i]))
    
    r  = requests.get(url)
    data = r.text
    soup = BeautifulSoup(data,features="html.parser")

    for link in soup.find_all('a'):
        if (filef == 1 and link.get('href').endswith(".csv")) or (filef == 2 and link.get('href').endswith(".root")) or (filef == 3 and (link.get('href').endswith(".csv") or link.get('href').endswith(".root"))):
            #print "Gets ", link.get('href')
            try:
                wget.download(url+"/"+link.get('href'),'/storage/eirikgr/POLA-03/')
            except:
                print("Could not find any data for " )
