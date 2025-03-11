__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

query_file = sys.argv[1]      # ITS1_seq_name_file.txt
result_file = sys.argv[2]	  # RES_ITS1.txt
sample_metadata = sys.argv[3] # GF5_SAMPLES_METADATA_TOTALFINAL_SAMPLE_ID.txt.gz
out_tab = sys.argv[4]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# CREATE SEQUENCE FILES
#############################################

metadata = {}
for line in openfile(sample_metadata):
     vals = line.rstrip().split('\t')
     metadata[vals[0]] = line.rstrip()
     #print(vals[0]+" "+vals[1])
print("Samples metadata loaded - " + str(len(metadata)))


queries = {}
for line in openfile(query_file):
     vals = line.rstrip().split('\t')
     queries[vals[0]] = vals[1]
     #print(vals[0]+" "+vals[1])
print("Queries loaded - " + str(len(queries)))

i = 0
fp = open(out_tab, "w")
out_line = "QUERY\tABUNDANCE\tsequence\tDatabase_ID\tPermanentID\tadd_date\tpaper_id\ttitle\tyear\tauthors\tjournal\tdoi\tcontact\tsample_name\tsample_type\tmanipulated\tsample_description\tlatitude\tlongitude\tcontinent\tyear_of_sampling\tBiome\tsequencing_platform\ttarget_gene\tprimers\tprimers_sequence\tsample_seqid\tsample_barcode\televation\tMAT\tMAP\tMAT_study\tMAP_study\tBiome_detail\tcountry\tmonth_of_sampling\tday_of_sampling\tplants_dominant\tplants_all\tarea_sampled\tnumber_of_subsamples\tsampling_info\tsample_depth\textraction_DNA_mass\textraction_DNA_size\textraction_DNA_method\ttotal_C_content\ttotal_N_content\torganic_matter_content\tpH\tpH_method\ttotal_Ca\ttotal_P\ttotal_K\tsample_info\tlocation\tarea_GPS\tITS1_extracted\tITS2_extracted\tITS_total\tpH_0cm_depth\tSOC_0cm_depth % \telev - GPS visualizer"
fp.write(out_line+"\n")
for line in openfile(result_file):
        l = line.rstrip()
        vals = l.split('\t')
        seq = vals[5]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        #
        if queries.has_key(seq):
            query = queries[seq]
        else:
            query = "[NOT EXACT]"
        #
        print("samples: "+str(len(sn)))
        for x in range(len(sn)):
            sample = sn[x] # get sample id
            abund = sa[x]
            fp.write(query + "\t" + str(abund) + "\t" + seq + "\t" + metadata[sample]+"\n")
            i += 1
fp.close()


print("DONE :) Records written: "+str(i))
