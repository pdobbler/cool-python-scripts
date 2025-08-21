__author__ = 'Wietrack 2023'

import sys
import operator
import hashlib
import gzip

samples_metadata = sys.argv[1]
samples_pairing = sys.argv[2]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

###########################
# Load pairing table
# 40590 N8833
###########################
sample_ids = {}
for line in openfile(samples_pairing):
    vals = line.strip().split("\t")
    sample_ids[vals[1]] = vals[0]


########################
# CHECKS
########################

platforms = ['DNBSEQ-G400', 'PacBio', 'Illumina', 'IonTorrent', '454Roche', 'Oxford Nanopore']

def check_platform(val):
    val = val.strip();
    if val in platforms:
        return val
    else:
        print("platform error: "+val)
        print(" "+0)
        return 0

markers = ['V4-V5', 'V4', 'V3-V4', 'V3-V5', '16S']

def check_marker(val):
    val = val.strip();
    if val in markers:
        return val
    else:
        print("marker error: "+val)
        print(" "+0)
        return 0

continents = ['Africa', 'Antarctica', 'Arctic Ocean', 'Asia', 'Atlantic Ocean', 'Australia', 'Europe', 'Indian Ocean', 'North America', 'Pacific Ocean', 'South America']

def check_continent(val):
    val = val.strip();
    if val in continents:
        return val
    else:
        print("continent error: "+val)
        print(" "+0)
        return 0

methods = ['KCl', 'NA_', 'H2O', 'CaCl2', 'in situ', 'KCl + CaCl2']

def check_ph_method(val):
    val = val.strip();
    if val in methods:
        return val
    else:
        print("method error: "+val)
        print(" "+0)
        return 0

def float_or_null(val):
    if val == "NA_":
        return "NULL"
    else:
        return str(float(val))

def int_or_null(val):
    if val == "NA_":
        return "NULL"
    else:
        return str(int(val))

def handle_int_range_with_null(val):
    if val == "NA_":
        return "NULL" + "\t" + "NULL"
    else:
        if "to" in val:
            #print(val.split("to")[0].strip() + "\t" +  val.split("to")[1].strip())
            return str(int(val.split("to")[0].strip())) + "\t" +  str(int(val.split("to")[1].strip()))
        else:
            return str(int(val)) + "\t" + str(int(val))

def handle_float_range_with_null(val):
    if val == "NA_":
        return "NULL" + "\t" + "NULL"
    else:
        if "to" in val:
            #print(val.split("to")[0].strip() + "\t" +  val.split("to")[1].strip())
            return str(float(val.split("to")[0].strip())) + "\t" +  str(float(val.split("to")[1].strip()))
        else:
            return str(float(val)) + "\t" + str(float(val))

def handle_depth_range_with_null(val):
    special = False
    if "/" in val:
        print("warning strange depth def: " + val)
        val = val.split("/")[0].split("to")[0].strip() + " to " + val.split("/")[1].split("to")[1].strip()
        special = True
    if "-" in val:
        print("warning strange depth def: " + val)
        val = val.split("-")[0] + " to " + val.split("-")[1] 
        special = True
    if val == "NA_":
        return "NULL" + "\t" + "NULL"
    else:
        if "to" in val:
            outval = str(float(val.split("to")[0].strip())) + "\t" +  str(float(val.split("to")[len(val.split("to")) - 1].strip()))
            if special:
                print("Corrected >>> "+outval)
            return outval
        else:
            return str(float(val)) + "\t" + str(float(val))

def handle_dna_mass_range_with_null(val):
    special = False
    #if "/" in val:
    #    print("warning strange depth def: " + val)
    #    val = val.split("/")[0].split("to")[0].strip() + " to " + val.split("/")[1].split("to")[1].strip()
    #    special = True
    #if "-" in val:
    #    print("warning strange depth def: " + val)
    #    val = val.split("-")[0] + " to " + val.split("-")[1] 
    #    special = True
    if val == "NA_":
        return "NULL" + "\t" + "NULL"
    else:
        if "to" in val:
            outval = str(float(val.split("to")[0].strip())) + "\t" +  str(float(val.split("to")[len(val.split("to")) - 1].strip()))
            if special:
                print("Corrected >>> "+outval)
            return outval
        else:
            return str(float(val)) + "\t" + str(float(val))


def handle_elevation(val,permanent_id):
    if val == "NA_":
        return "NULL"
    else:
        if "." in val:
            #print("Warning: Sample " + permanent_id + " - elevation : " + val + " -> " + str(int(round(float(val)))))
            return str(int(round(float(val))))
        else:
            return str(int(val))

#def handle_area_GPS(val):
#    if "to" in val:
#        return val.split("to")[0].strip() + "-" +  val.split("to")[1].strip()
#    else:
#        if val == "NA_":
#            return "NULL"
#        else:
#            return str(float(val))

def handle_chemism(val):
    if "<" in val:
        val = val[1:]
    if val == "NA_":
        return "NULL"
    else:
        return str(float(val))

def boolean_to_tinyint(val):
    val = val.lower()
    if val == "false":
        return 0
    if val == "true":
        return 1
    print("Error boolean "+val)
    print(" "+0)



##################################
# permanent_id,[0]
#`add_date` varchar(10) NOT NULL,[1]            :P (advanced)
#`paper_id`varchar(32)NOTNULL,[2]               :)
#`title` TEXT NOT NULL,[3]                      :)
#`year` varchar(4) NOT NULL,[4]                 :)
#  `authors` TEXT NOT NULL,[5]                  :)
#  `journal` TEXT NOT NULL,[6]                  :)
#  `doi` TEXT NOT NULL,[7]                      :)
#  `contact` TEXT NOT NULL,[8]                  :)
#  `sample_name` TEXT NOT NULL,[9]              :P (advanced)
#  `sample_type` TEXT NOT NULL,[10]             ;D
#`manipulated` varchar(5) NOT NULL,[11]         :)
#  `sample_description` TEXT NOT NULL,[12]      :P (advanced)
#  `latitude` float NOT NULL,[13]               :)  -> must be float!! Not empty!
#  `longitude` float NOT NULL,[14]              :)  -> must be float!! Not empty!
#  `continent` varchar(32) NOT NULL,[15]        :P (advanced)
#  `year_of_sampling` varchar(32) NOT NULL,[16] :)  -> must be int!! Not empty!
#  `Biome` TEXT NOT NULL,[17]                   :)
#`sequencing_platform`varchar(32)NOTNULL,[18]   :P (advanced)
#`target_gene`varchar(32)NOTNULL,[19]           :P (advanced)
#  `primers` TEXT NOT NULL,[20]                 :)
#  `primers_sequence` TEXT NOT NULL,[21]        :P (advanced)
#  `sample_seqid` TEXT NOT NULL,[22]            :P (advanced)
#  `sample_barcode` TEXT NOT NULL,[23]          :P (advanced)
#  `elevation` varchar(32) NOT NULL,[24]        :P (advanced) -> must be int!! Not empty!
#  `MAT` varchar(32) NOT NULL,[25]              :)  -> must be float!! Not empty!
#  `MAP` varchar(32) NOT NULL,[26]              :)  -> must be float!! Not empty!
#  `MAT_study` varchar(32) NOT NULL,[27]        :P (advanced)
#  `MAP_study` varchar(32) NOT NULL,[28]        :P (advanced)
#  `Biome_detail` TEXT NOT NULL,[29]            :P (advanced)
#  `country` TEXT NOT NULL,[30]                 :P (advanced)
#`month_of_sampling`varchar(32)NOTNULL,[31]     :P (advanced)
#  `day_of_sampling` varchar(32) NOT NULL,[32]  :P (advanced)
#  `plants_dominant` TEXT NOT NULL,[33]         :P (advanced)
#  `plants_all` TEXT NOT NULL,[34]              :P (advanced)
#  `area_sampled` varchar(32) NOT NULL,[35]     :P (advanced)
#  `number_of_subsamples` varchar(32) NOT NULL,[36] :P (advanced)
#  `sampling_info` TEXT NOT NULL,[37]           :P (advanced)
#  `sample_depth` varchar(32) NOT NULL,[38]     :P (advanced)
#  `extraction_DNA_mass` varchar(32) NOT NULL,[39] :P (advanced)
#  `extraction_DNA_size` TEXT NOT NULL,[40]     :P (advanced)
#  `extraction_DNA_method` TEXT NOT NULL,[41]   :P (advanced)
#  `total_C_content` varchar(32) NOT NULL,[42]  :P (advanced)
#  `total_N_content` varchar(32) NOT NULL,[43]  :P (advanced)
#  `organic_matter_content` varchar(32) NOT NULL,[44] :P (advanced)
#  `pH` varchar(32) NOT NULL,[45]               :)
#  `pH_method` varchar(64) NOT NULL,[46]        :P (advanced)
#  `total_Ca` varchar(32) NOT NULL,[47]         :P (advanced)
#  `total_P` varchar(32) NOT NULL,[48]          :P (advanced)
#  `total_K` varchar(32) NOT NULL,[49]          :P (advanced)
#  `sample_info` TEXT NOT NULL,[50]             :P (advanced)
#  `location` TEXT NOT NULL,[51]                :P (advanced)
#  `area_GPS` varchar(32) NOT NULL,[52]         :P (advanced)
#  `ITS1_extracted` int NOT NULL,[53]           :P (advanced)
#  `ITS2_extracted` int NOT NULL,[54]           :P (advanced)
#  `ITS_total` int NOT NULL[55]                 :)
#  pH_0cm_depth [56]  
#  SOC_0cm_depth % [57]  
#  elev - GPS visualizer[58]  
##################################

study_id = 0
current_study = 0
studies_check = {}
studies = []
samples_basic = {}
samples_advanced = {}
coords = {}
paper_check = {}
paper_manipulated = {}

#studies["tadd_date\tpaper_id\tyear\ttitle\tauthors\tjournal\tdoi\tcontact\tmanipulated"] = 0 
#samples_basic[0] = "permanent_id\tsample_type\tlatitude\tlongitude\tyear_of_sampling_from\tyear_of_sampling_to\tbiome\tprimers\tMAT\tMAP\tpH\tITS_total\tmanipulated" 
#samples_advanced[0] = ("paper\tadd_date\tsample_name\tsample_description\tcontinent\tsequencing_platform\ttarget_gene\tprimers_sequence\tsample_seqid\tsample_barcode\televation\tMAT_study\tMAP_study\t" 
#+ "Biome_detail\tcountry\tmonth_of_sampling\tday_of_sampling\tplants_dominant\tplants_all\tarea_sampled\tnumber_of_subsamples_from\tnumber_of_subsamples_to\tsampling_info\tsample_depth_from\tsample_depth_to\textraction_DNA_mass_from\textraction_DNA_mass_to\textraction_DNA_size\t" 
#+ "extraction_DNA_method\ttotal_C_content\ttotal_N_content\torganic_matter_content\tpH_study\tpH_method\ttotal_Ca\ttotal_P\ttotal_K\tsample_info\tlocation\tarea_GPS_from\tarea_GPS_to\tITS1_extracted\tITS2_extracted")
i = 1
for line in openfile(samples_metadata):
    vals = line.strip().split("\t")
    # studies
    if paper_check.has_key(vals[2]):
        if vals[1] != paper_check[vals[2]]:
            vals[1] = paper_check[vals[2]]
    else:
        paper_check[vals[2]] = vals[1]
    #
    study_str =   vals[1] + "\t" + vals[2] + "\t" + str(int(vals[4])) + "\t" + vals[3] + "\t" + vals[5] + "\t" + vals[6] + "\t" + vals[7] + "\t" + vals[8]
    if not studies_check.has_key(study_str):
        study_id += 1
        studies_check[study_str] = study_id
        studies.append(str(study_id) + "\t" + study_str)
        current_study = study_id
        paper_manipulated[study_id] = boolean_to_tinyint(vals[11])
    else:
        # check manupulated        
        current_study = studies_check[study_str]
        if boolean_to_tinyint(vals[11]) == 1:
            paper_manipulated[current_study] = 1

    # samples 
    print(vals[0]+" "+vals[38]+" "+str(i))
    samples_basic[sample_ids[vals[0]]] = (str(current_study) + "\t" + vals[0] +  "\t" + vals[10] + "\t" + str(float(vals[13])) + "\t" + str(float(vals[14])) + "\t" + check_continent(vals[15]) + "\t" + handle_int_range_with_null(vals[16]) 
    + "\t" + vals[17] + "\t" + vals[20] 
    + "\t" + float_or_null(vals[25]) 
    + "\t" + float_or_null(vals[26]) 
    + "\t" + float_or_null(vals[56])  
    + "\t" + float_or_null(vals[57]) 
    + "\t" + str(int(vals[55])) 
    + "\t" + str(boolean_to_tinyint(vals[11])))
    #print(samples_basic[sample_ids[vals[0]]])
    # samples advanced
    samples_advanced[sample_ids[vals[0]]] = (vals[9] + "\t" + vals[12] + "\t" + check_platform(vals[18]) + "\t" + check_marker(vals[19]) + "\t" + vals[21] + "\t" + vals[22] 
    + "\t" + vals[23] + "\t" + handle_elevation(vals[58],vals[0]) + "\t" + float_or_null(vals[27]) + "\t" + float_or_null(vals[28]) + "\t" + vals[29] + "\t" + vals[30] + "\t" + vals[31] + "\t" + vals[32] + "\t" + vals[33]
    + "\t" + vals[34] + "\t" + float_or_null(vals[35]) + "\t" + handle_int_range_with_null(vals[36]) + "\t" + vals[37] + "\t" + handle_depth_range_with_null(vals[38]) + "\t" + handle_dna_mass_range_with_null(vals[39]) + "\t" + vals[40] + "\t" + vals[41] + "\t" 
    + handle_chemism(vals[42]) + "\t" + float_or_null(vals[43]) + "\t" + float_or_null(vals[44]) + "\t" + float_or_null(vals[45]) + "\t" + check_ph_method(vals[46]) + "\t" + float_or_null(vals[47]) + "\t" + handle_chemism(vals[48]) + "\t" + float_or_null(vals[49]) 
    + "\t" + vals[50] 
    + "\t" + vals[51] 
    + "\t" + handle_float_range_with_null(vals[52]) 
    + "\t" + str(int(vals[53])) 
    + "\t" + str(int(vals[54])))
    i += 1

print("Studies: "+str(len(studies)))

# write tables

#################################################################################

pf = open("SAMPLES_PAPERS.txt", 'w')

# Sort the `studies` dictionary by values in ascending order
for s in studies:
    vals = s.split("\t")
    #print(str(sorted_studies[study_str]))
    pf.write(vals[0] + "\t" + vals[1]  + "\t" + vals[3]  + "\t" + vals[4]  + "\t" + vals[5] + "\t" + vals[6] + "\t" + vals[7]  + "\t" + vals[8] + "\t" + str(paper_manipulated[int(vals[0])]) + "\n")
pf.close()

#  `id` int NOT NULL,
#  `add_date` VARCHAR(10) NOT NULL,
#  `paper_id` VARCHAR(48) NOT NULL,
#  `year` int NOT NULL,
#  `title` TEXT NOT NULL,
#  `authors` TEXT NOT NULL,
#  `journal` VARCHAR(128) NOT NULL,
#  `doi` VARCHAR(64) NOT NULL,
#  `contact` TEXT NOT NULL,
#  `manipulated` TINYINT(1) NOT NULL

pf = open("SAMPLES_PAPERS_SQL.txt", 'w')
i = 0
header = ("INSERT INTO `samples_papers` (`id`, `add_date`, `year`, `title`, `authors`, `journal`,  `doi`, `contact`,  `manipulated`) VALUES "+ "\n")
pf.write(header)
for s in studies:
    i += 1
    vals = s.split("\t")
    text = ("(" + vals[0] + ",'" + vals[1] + "'," + vals[3]  + ",'" 
        + vals[4]  + "','" + vals[5] + "','" + vals[6] + "','" + vals[7]  + "','" + vals[8] +"'," + str(paper_manipulated[int(vals[0])])  + ")")
    if i == len(studies):
        pf.write(text+";\n")
    else:
        pf.write(text+",\n")
pf.close()

#################################################################################

pf = open("SAMPLES_BASIC.txt", 'w')
for key in samples_basic:
    #print(str(key))
    pf.write(str(key) + "\t" + samples_basic[key] + "\n")
pf.close()

#CREATE TABLE IF NOT EXISTS `samples_basic` (
#  `id` int NOT NULL,
#  `permanent_id` VARCHAR(11) NOT NULL,
#  `sample_type` VARCHAR(32) NOT NULL,
#  `latitude` float NOT NULL,
#  `longitude` float NOT NULL,
#  `year_of_sampling_from` int NOT NULL,
#  `year_of_sampling_to` int NOT NULL,
#  `Biome` VARCHAR(32) NOT NULL,
#  `primers` VARCHAR(128) NOT NULL,
#  `MAT` FLOAT,
#  `MAP` FLOAT,
#  `pH` FLOAT,
#  `ITS_total` int NOT NULL,
#  `manipulated` TINYINT(1) NOT NULL
#);

pf = open("SAMPLES_BASIC_SQL.txt", 'w')
count = 0
i = 0
header = ("INSERT INTO `samples_basic` (`id`, `paper`, `permanent_id`, `sample_type`, `latitude`, `longitude`, `continent`, `year_of_sampling_from`,"
  +"`year_of_sampling_to`,  `Biome`, `primers`,  `MAT`, `MAP`,  `pH`, `SOC`, `ITS_total`,  `manipulated`) VALUES "+ "\n")
pf.write(header)
for key in samples_basic:
    count += 1
    i += 1
    vals = samples_basic[key].split("\t")
    text = ("(" + str(key) + ","+vals[0]+ ",'"+vals[1]+ "','"+vals[2]+ "',"+vals[3]+ ","+vals[4] + ",'" + vals[5]+ "',"+vals[6]+ 
     ","+vals[7]+ ",'" +vals[8]+ "','"+vals[9]+ "',"+vals[10]+ "," +vals[11]+ "," +vals[12]+ "," +vals[13]+ "," +vals[14] + "," + vals[15]  + ")")
    if i == len(samples_basic):
        pf.write(text+";\n")
    else:
        if count>10000:
            count = 0
            pf.write(text+";\n")
            pf.write(header)
        else:
            pf.write(text+",\n")
pf.close()

#################################################################################

pf = open("SAMPLES_ADVANCED.txt", 'w')
for key in samples_advanced:
    #print(str(key))
    pf.write(str(key) + "\t" + samples_advanced[key] + "\n")
pf.close()


pf = open("SAMPLES_ADVANCED_SQL.txt", 'w')
count = 0
i = 0
header = ("INSERT INTO `samples_advanced` (`id`, `sample_name`, `sample_description`, `sequencing_platform`,"
  +"`target_gene`,  `primers_sequence`, `sample_seqid`,  `sample_barcode`,"
  +"`elevation`,  `MAT_study`, `MAP_study`,  `Biome_detail`, `country`,  `month_of_sampling`, `day_of_sampling`, `plants_dominant`,"
  +"`plants_all`, `area_sampled`, `number_of_subsamples_from`, `number_of_subsamples_to`, `sampling_info`, `sample_depth_from`,  `sample_depth_to`,"
  +"`extraction_DNA_mass_from`, `extraction_DNA_mass_to`, `extraction_DNA_size`, `extraction_DNA_method`, `total_C_content`, `total_N_content`,"
  +"`organic_matter_content`, `pH_study`, `pH_method`, `total_Ca`, `total_P`, `total_K`, `sample_info`, `location`, `area_GPS_from`, `area_GPS_to`, `ITS1_extracted`, `ITS2_extracted`"
+") VALUES "+ "\n")
pf.write(header)
for key in samples_advanced:
    #print(str(key))
    count += 1
    i += 1
    vals = samples_advanced[key].split("\t")
    text = ("(" + str(key) + ",'" + vals[0] + "','"+vals[1]+ "','"+vals[2]+ "','"+vals[3]+ "','"+vals[4]+ "','" + vals[5]+ "','"+ vals[6]+ "'," 
        + vals[7]+"," + vals[8]+"," + vals[9] 
        + ",'" + vals[10] + "','" + vals[11] + "','"+vals[12] + "','"+vals[13] + "','"+vals[14] + "','"+vals[15] + "',"
        +vals[16] + ","+vals[17]+ "," + vals[18] 
        + ",'" + vals[19] + "'," 
        +vals[20] + "," +vals[21] + "," +vals[22] + "," +vals[23] 
        + ",'" +vals[24] + "','" +vals[25] + "'," 
        +vals[26] + "," +vals[27] + ","  + vals[28] + "," + vals[29] 
        + ",'" +vals[30] + "'," 
        +vals[31] + "," +vals[32] + ","+vals[33] 
        + ",'"+vals[34] + "','"+vals[35] + "',"
        +vals[36] + ","+vals[37] + ","+vals[38] + "," + vals[39] 
        + ")")
    if i == len(samples_advanced):
        pf.write(text+";\n")
    else:
        if count>10000:
            count = 0
            pf.write(text+";\n")
            pf.write(header)
        else:
            pf.write(text+",\n")
pf.close()

print("Done :)")







