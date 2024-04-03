__author__ = 'vetrot'

import sys
import urllib
import urllib2
import re

input_file = sys.argv[1]
input_type = sys.argv[2] #A - accessions, N - names
taxa_output = sys.argv[3]


if input_type == 'A':
    print('get taxonomy by ACCESSION numbers!')
else:
    print('get taxonomy by TAXON names!')

query = []
taxa_unique = {}
#process file...
for line in open(input_file):
    val = line.rstrip()
    query.append(val)
    taxa_unique[val] = {}

print(str(len(query))+' queries loaded - resulting in '+str(len(taxa_unique))+' unique queries.')
#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=Russula[Name]All[Organism]&retmode=xml

def save_taxonomy(fp, taxa_unique, name):
    taxa_ranks = taxa_unique[name]
    line = 'no_hit' + '\t' + 'no_hit' + '\t' + 'no_hit' + '\t' + 'no_hit' + '\t' + 'no_hit' + '\t' + 'no_hit' + '\t' + 'no_hit' + '\t[' + name + ']'
    if taxa_ranks.has_key('superkingdom'):
        # superkingdom - domain
        superkingdomTax = taxa_ranks['superkingdom']
        # phylum
        phylumTax = 'undefined ' + superkingdomTax
        if taxa_ranks.has_key('phylum'):
            phylumTax = taxa_ranks['phylum']
        # class
        classTax = 'undefined ' + phylumTax.replace('undefined ', '')
        if taxa_ranks.has_key('class'):
            classTax = taxa_ranks['class']
        # order
        orderTax = 'undefined ' + classTax.replace('undefined ', '')
        if taxa_ranks.has_key('order'):
            orderTax = taxa_ranks['order']
        # family
        familyTax = 'undefined ' + orderTax.replace('undefined ', '')
        if taxa_ranks.has_key('family'):
            familyTax = taxa_ranks['family']
        # family
        genusTax = 'undefined ' + familyTax.replace('undefined ', '')
        if taxa_ranks.has_key('genus'):
            genusTax = taxa_ranks['genus']
        # organism
        organism = name
        if taxa_ranks.has_key('species'):
            organism = taxa_ranks['species']
        line = superkingdomTax + '\t' + phylumTax + '\t' + classTax + '\t' + orderTax + '\t' + familyTax + '\t' + genusTax + '\t' + organism + '\t[' + name + ']'
    #print(line)
    fp.write(line + '\n')

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def find_between_r( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""

fp = open('taxa_processed.txt', 'w')
for taxName in taxa_unique:
    if input_type == 'N':
        # taxon names
        url_taxa = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term='+taxName+'[Name]All[Organism]&retmode=xml'
        url_taxa = url_taxa.replace(" ", "%20")
        #print (url_taxa)
        #file = urllib2.urlopen(url_taxa)
        file = urllib.urlopen(url_taxa)
        data = file.read()
        file.close()
        taxid = find_between(data, "<Id>", "</Id>")
    else:
        # try nucleotide
        url_taxa = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + taxName + '&rettype=fasta&retmode=xml'
        file = urllib.urlopen(url_taxa)
        data = file.read()
        file.close()
        taxid = find_between(data, "<TSeq_taxid>", "</TSeq_taxid>")
        if len(taxid) == 0:
            # try protein
            url_taxa = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=' + taxName + '&rettype=fasta&retmode=xml'
            file = urllib.urlopen(url_taxa)
            data = file.read()
            file.close()
            taxid = find_between(data, "<TSeq_taxid>", "</TSeq_taxid>")

    print(taxName+' > TAXID >'+taxid)
    if len(taxid) > 0:
        data = ''
        iter = 0
        while data == '':
            print('iteration: '+str(iter))
            url_taxid = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id='+taxid+'&retmode=xml'
            #file = urllib2.urlopen(url_taxid)
            file = urllib.urlopen(url_taxid)
            data = file.read()
            file.close()
            iter += 1

        #print(data)
        taxa_ranks = {}
        parts = data.split('<Taxon>')
        if len(parts) > 0:
            for i in range(0,len(parts)):
                #print(parts[i])
                rank = find_between(parts[i], "<Rank>", "</Rank>")
                #print (rank)
                taxa_ranks[rank] = find_between(parts[i], "<ScientificName>", "</ScientificName>")
        taxa_unique[taxName] = taxa_ranks
        save_taxonomy(fp, taxa_unique, taxName)
    #break

fp.close()

print ('Taxonomy download finished...')

fp = open(taxa_output, 'w')
line = 'domain\tphylum\tclass\torder\tfamily\tgenus\torganism\ttax_key'
#print(line)
fp.write(line + '\n')
for name in query:
    save_taxonomy(fp, taxa_unique, name)

fp.close()
print('Done :)')




