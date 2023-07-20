from Bio import SeqIO, Entrez
import pandas as pd

#searchTerm = '(("colorectal cancer") OR ("rectal cancer") OR CRC OR ("bowel cancer") OR ("colon cancer") OR adenoma) AND ( (shotgun OR ("shotgun metagenomics") OR WGS OR ("metagenomics") OR ("metagenomic") OR ("shotgun sequencing") OR ("shotgun metagenomics sequencing")) OR (16S OR amplicon OR ("ribosomal sequencing") ))'
#searchTerm = '(("colorectal cancer") OR ("rectal cancer") OR CRC OR ("bowel cancer") OR ("colon cancer") OR adenoma) AND ((shotgun OR ("shotgun metagenomics") OR WGS OR ("metagenomics") OR ("metagenomic") OR ("shotgun sequencing") OR ("shotgun 	metagenomics sequencing")) OR ( (16S OR amplicon OR "ribosomal" OR "ribosome" ) AND "sequencing") OR ( ("mucosa" or "mucosal" or "mucosae" or "tissue") AND ("sequencing") ))'
#searchTerm = '(("colorectal cancer") OR ("rectal cancer") OR CRC OR ("bowel cancer") OR ("colon cancer") OR adenoma) AND (microbiome OR microbiota OR community OR bacteria OR bacterial OR microbial) ( ( (metagenome OR metagenomes OR "metagenome-wide" OR shotgun OR ("shotgun metagenomics") OR WGS OR ("metagenomics") OR ("metagenomic") OR ("shotgun") AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile") ) OR  ( (16S OR amplicon OR "ribosomal" OR "ribosome" ) AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile")) OR ( ("mucosa" or "mucosal" or "mucosae" or "tissue") AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile") ) )'
#searchTerm = '(("colorectal cancer") OR ("rectal cancer") OR CRC OR ("bowel cancer") OR ("colon cancer") OR adenoma) ( ( (metagenome OR metagenomes OR "metagenome-wide" OR shotgun OR ("shotgun metagenomics") OR WGS OR ("metagenomics") OR ("metagenomic") OR ("shotgun") AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile" OR "profiling") ) OR  ( (16S OR amplicon OR "ribosomal" OR "ribosome" ) AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile" OR "profiling")) OR ( ("mucosa" or "mucosal" or "mucosae" or "tissue") AND ("sequencing" OR "sequenced" OR "catalog" OR "catalogue" OR "catalogued" OR "profiled" OR "profile" OR "profiling") ) )'
#searchTerm = '(("colorectal cancer") OR ("rectal cancer" OR CRC OR ("bowel cancer") OR ("colon cancer") OR adenoma) AND ( (shotgun OR metagenome OR metagenomes OR metagenomics OR (meta AND sequencing) OR (meta-sequencing) OR ("shotgun metagenomics") OR WGS OR ("metagenomics") OR ("metagenomic") OR ("shotgun sequencing") OR ("shotgun metagenomics sequencing")) OR (16S OR amplicon OR ("ribosomal sequencing") ))'
#searchTerm = '( ( (colon OR colorectal or rectal or bowel or intest*) AND (cancer OR polyp OR adenoma OR neoplas* OR malign*) ) AND ( (shotgun OR metagenom* OR (meta AND sequencing) OR WGS) OR (16S OR amplicon OR ("ribosomal sequencing") )  OR (microb* AND signature*) OR (microb* AND profil*) OR ("associated microb"*) )'
#searchTerm='( (colon OR colorectal or rectal or bowel) AND (cancer OR polyp OR adenoma OR neoplas* OR malign*) ) AND ( ( shotgun OR metagenom* OR (meta AND sequencing) OR WGS ) OR (16S OR amplicon OR ("ribosomal sequencing") ) OR (microb* AND signature*) OR ("associated microb*") OR ("*associated bacteri*") OR (profile* AND microbi*) OR "bacterial taxon markers")'
#searchTerm = '( (colon OR colorectal or rectal or bowel or intest*) AND (cancer OR polyp OR adenoma OR neoplas* OR malign*) ) AND ( shotgun OR metagenom* OR (meta AND sequencing) OR WGS OR 16S OR amplicon OR "ribosomal sequencing" OR (microb* AND signature*) OR "associated microb*" OR "*associated bacteri*" OR (profile* AND microbi*) OR "bacterial taxon markers" ) AND ( (human* OR patient* OR participant*) NOT (mouse OR mice OR murin* OR rat OR rats) )'
#searchTerm = '( (colon[tiab] OR colorectal[tiab] or rectal[tiab] or bowel[tiab] or intest*[tiab]) AND (cancer[tiab] OR carcinoma[tiab] OR polyp[tiab] OR adenoma[tiab] OR neoplas*[tiab] OR malign*[tiab]) ) AND ( shotgun[tiab] OR metagenom*[tiab] OR (meta[tiab] AND sequencing[tiab]) OR WGS[tiab] OR 16S[tiab] OR amplicon[tiab] OR "ribosomal sequencing"[tiab] OR (microb*[tiab] AND signature*[tiab]) OR "associated microb*"[tiab] OR "*associated bacteri*"[tiab] OR (profile*[tiab] AND microbi*[tiab]) OR "bacterial taxon markers"[tiab] OR "faecal bacterial taxa"[tiab] OR "Analysis of the gut microbiome from stool samples" ) AND ( (human*[tiab] OR patient*[tiab] OR participant*[tiab] OR subject*[tiab] OR cohort*[tiab] OR individual*[tiab] OR "tumor/normal pair*"[tiab]) )'
searchTerm = '( (colon[tiab] OR colorectal[tiab] OR rectal[tiab] OR bowel[tiab]) AND (cancer*[tiab] OR precancer*[tiab] OR carcino*[tiab] OR polyp[tiab] OR adenoma*[tiab] OR neoplas*[tiab] OR malign*[tiab]) ) AND ( shotgun[tiab] OR metagenom*[tiab] OR (meta[tiab] AND sequencing[tiab]) OR WGS[tiab] OR 16S[tiab] OR amplicon[tiab] OR "ribosomal sequencing"[tiab] OR "associated microb*"[tiab] OR "*associated bacteri*"[tiab] OR (profile*[tiab] AND microbi*[tiab]) OR "bacterial tax*"[tiab] ) AND ( (human*[tiab] OR patient*[tiab] OR participant*[tiab] OR subject*[tiab] OR cohort*[tiab] OR individual*[tiab] OR "tumor/normal pair*"[tiab]) ) NOT ( "Case Reports" [publication type]) NOT ( "Review"[publication type] ) AND (2011:2023[dp])'
Entrez.email = 'nicolai.karcher@embl.de'

res = Entrez.esearch(db = 'PubMed', term = searchTerm, retmax = 10000000)
res = Entrez.read(res)

print(len(res['IdList']))
a = len(res['IdList'])
#records = Entrez.esummary(db="PubMed", id=",".join(res['IdList']), report="full")
for i in range(5):
    print("Try number {}".format(i)) 
    try:
        records = Entrez.efetch(db="PubMed", id=",".join(res['IdList']), report="full")
        records = Entrez.read(records)
        break
    except Exception as e:
        print(e)
        continue

print('...here')
# Most papers are returned in the 'PubmedArticle' field, very very few are returned 'PubmedBookArticle'.
b = len(records['PubmedArticle'])

if a != b:
    print("I cant find all articles in pubmed")
    print("if the following numbers are very similar, dont mind")
    print(a)
    print(b)

dat = []
for record in records['PubmedArticle']:
    d = record['MedlineCitation']['Article']
    try:
        pmid = record['MedlineCitation']["PMID"]
    except:
        pmid = None
    try:
        date = d['ArticleDate']
    except:
        date = None
    try:
        date2 = d['Journal']['JournalIssue']['PubDate']['Year']
    except:
        date2 = None
    try:
        authors = d['AuthorList']
        firstAuthor = d['AuthorList'][0]['LastName']
    except:
        authors = None
        firstAuthor
    try:
        title = d['ArticleTitle']
    except:
        title = None
    try:
        abstract = " ".join([str(x) for x in d['Abstract']['AbstractText']])
    except:
        abstract = None
    try:
        pubType = d['PublicationTypeList']
    except:
        pubType = None        
    dat.append([date, date2, firstAuthor, title, abstract, pubType, pmid])

dataParsed = pd.DataFrame(dat)
dataParsed = dataParsed.rename(columns={0: "date", 1: "dateAlternative", 2: "firstAuthor", 3 : 'title', 4 : "abstract", 5:'pubType', 6:'PMID'})

# IMPORTANT: The pubmed query searches for strings in double hyphens (such as "bacterial taxon markers") as exact string matches. However, this doesn't always work, but you only see that it doesn't work when you
# Manually search in the pubmed web GUI.


# a few studies have no abstract, so exclude
#dataParsed['hasAbstract'] = [True if x is not None else False for x in dataParsed['abstract']]
#dataParsed = dataParsed[dataParsed['hasAbstract']]

# Exclude potential reviews 
dataParsed['potentialReview'] = [True if ("review" in x or "Review" in x) else False for x in dataParsed['pubType']] 
#dataParsed = dataParsed[[False if x else True for x in dataParsed['potentialReview']]]

# Consider only stuff newer than 2011. dateAlternative is more reliable and less sparse than date.
#dataParsed = dataParsed[[True if int(x) >= 2012 else False for x in dataParsed['dateAlternative']]]

# Mark studies with 'mice' in abstract
dataParsed['mice'] = [True if ('mice' in x or "Mice" in x or "murine" in x or "murin" in x or "Murine" in x or "Murin" in x) else False for x in dataParsed['abstract']]

# Mark studies with 'mice' in abstract
dataParsed['human'] = [True if ('human' in x or "Human" in x or "patient" in x or "Patient" in x or 'patients' in x or 'Patients' in x) else False for x in dataParsed['abstract']]
#dataParsed = dataParsed[dataParsed['human']]

dataParsed.to_csv('new_DSs.tsv', sep = "\t")
