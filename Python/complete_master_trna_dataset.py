

"""
The purpose of this file is to collect the SGD IDs for each of the tRNA genes from the GtRNAdb, since for some reason
they are not included in the main table from that website. Since I'm at it, I will collect the RNAcentral IDs (even
though I don't know what they are) just in case. Then grab the systematic gene IDs and standard gene IDs from the SGD
database and add them to this as well. Then add all of this information to the master_dataset I generated in R before
(at "C:/MyStuff/tRNAs/Scripts/R/Mine/2.create_master_dataframe.R).
"""

# How about adding a little change to here, what happens with Git?
# This is another change, does this update automatically or what? Okay making 1 more change to see if now I can commit


#TODO
"""
Honestly I think it would be better to put this all together in a Python script, what I do in that R file and what I do 
here. I tried to do it quickly but the downloading the main gtrnadb dataset from Python didn't work :(. So yeah, looks 
like that won't be happening for now. 
"""



# 0. Set-up
## 0.1. Libraries
import urllib.request
from html_table_parser.parser import HTMLTableParser
import pandas as pd
import requests

## 0.2. Details to be able to access websites from Python
# Make sure I can connect to the Internet through the proxy
proxies = {"http": 'http://algo12:Adult_Garment72!@proxy.charite.de:8080',
           "https": 'http://algo12:Adult_Garment72!@proxy.charite.de:8080'}

# Add these headers just in case the website requires them, even though I don't think it does
headers = {'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:134.0) Gecko/20100101 Firefox/134.0",
           'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
           'Accept-Language': 'zh-CN,zh;q=0.9,en;q=0.8,ja;q=0.7'}

# Set up a function that collects the html content from a website (opens a website and read its binary contents (HTTP Response Body))
def url_get_contents(url):
    # making request to the website
    req = urllib.request.Request(url=url, headers=headers, timeout=60, verify=False)
    f = urllib.request.urlopen(req)

    # reading contents of the website
    return f.read()



# 1. Collect the tables
## 1.1. Collect the names of the tRNAs from the master dataset in order to build the URLs from them
gtrnas = pd.read_csv("C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\master_tRNA_dataset.csv")
gene_symbols = gtrnas["GtRNAdb_gene_symbol"].to_list()


## 1.2. Set up a function that collects the html content from a website
# Opens a website and read its binary contents (HTTP Response Body)
def url_get_contents(url):
    # making request to the website
    req = urllib.request.Request(url=url, headers=headers, timeout=60, verify=False)
    f = urllib.request.urlopen(req)

    # reading contents of the website
    return f.read()


## 1.3. Iterate over the tRNA genes, downloading the table for each and collecting its SGD ID
SGD_IDs = []
RNAcentral_IDs = []
for symbol in gene_symbols:
    # Build the URL
    url = "https://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/genes/" + symbol + ".html"

    # Get the html contents of the URL
    xhtml = requests.get(url, proxies=proxies, headers=headers, timeout=60, verify=False).content.decode('utf-8')

    # Defining the HTMLTableParser object
    p = HTMLTableParser()

    # feeding the html contents in the HTMLTableParser object
    p.feed(xhtml)

    # Now finally obtaining the data of the table required
    df = pd.DataFrame(p.tables[0], columns=["Parameter", "Value"])
    SGD_IDs.append(df.loc[df["Parameter"] == "SGD ID", "Value"].iloc[0])
    RNAcentral_IDs.append(df.loc[df["Parameter"] == "RNAcentral ID", "Value"].iloc[0])


## 1.4. Construct a table from these vectors, that will be the base for the output we produce here
data = [gene_symbols, SGD_IDs, RNAcentral_IDs]
out = pd.DataFrame(data).T
out.columns = ["GtRNAdb_gene_symbol", "Gene.primaryIdentifier", "RNAcentral_IDs"]



## 2. Load SGD table and add new columns - standard and systematic gene names, based on SGD ID
SGD = pd.read_csv("C:\\MyStuff\\ScRAP\\Code\\Biological questions\\Data\\yeastmine_results_2.tsv", sep="\t")
out = out.merge(SGD, on="Gene.primaryIdentifier", how="left")
out = out[["GtRNAdb_gene_symbol", "Gene.primaryIdentifier", "RNAcentral_IDs", "Gene.secondaryIdentifier", "Gene.symbol"]]



## 3. Add the columns in this dataset to the original master_dataset (created in R) and rewrite this one
master_dataset = pd.read_csv("C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\master_tRNA_dataset.csv")
master_dataset = master_dataset.merge(out, on="GtRNAdb_gene_symbol", how="left")
master_dataset.to_csv("C:\\MyStuff\\tRNAs\\Data\\GtRNAdb\\master_tRNA_dataset.csv", header=True, index=False)





























