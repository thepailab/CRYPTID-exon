#!/usr/bin/env python

import pandas as pd
import numpy as np
import ruptures as rpt
import os 
import argparse


### import arguments 
parser = argparse.ArgumentParser()
parser.add_argument('-splice_file', required=True,help = 'Output file from CRYPTID-SS or text file containing splice junctions (.csv).')
parser.add_argument('-bam', required=True,help = 'BAM file for exon prediction')
parser.add_argument('-intron_window', type=int, default= 100, help = 'size of coverage window extended into the intron from exon boundary. Default is 100 nt ')
parser.add_argument('-exon_window', type=int, default=200, help = 'size of coverage window extended into the exon from exon boundary. Default is 200 nt')
parser.add_argument('-output_dir', required=True, help = 'path to output directory')
args = parser.parse_args()


# parsing argument 
juncs_file=args.splice_file
bam=args.bam
intron_window = args.intron_window
exon_window = args.exon_window
outpath=args.output_dir

print("importing splice junctions...")
# import the splice junctions text file 
juncs = pd.read_csv(juncs_file, sep =",", low_memory = False)

# rename the needed columns 
juncs.rename(columns=dict(zip(juncs.columns[:5], ['jnc_id','chr','ss5','ss3','strand'])), inplace=True)


print("processing splice junctions...")

# make IDs for ss3 and ss5 in each junctions and make list of plus and minus ss3 and ss5.
juncs["id_ss3"] = juncs["chr"].astype(str) + "_" + juncs["ss3"].astype(str) + "_" + juncs["strand"].astype(str)
juncs["id_ss5"] = juncs["chr"].astype(str) + "_" + juncs["ss5"].astype(str) + "_" + juncs["strand"].astype(str)

# divide in plus and minus junctions 
juncs_plus = juncs[juncs["strand"] == "+"]
juncs_minus = juncs[juncs["strand"] == "-"]

# add exon start and exon end 
juncs["exon_start"] = "NA"
juncs.loc[juncs["strand"] == "+", "exon_start"] = juncs.loc[juncs["strand"] == "+", "ss3"] +1
juncs.loc[juncs["strand"] == "-", "exon_start"] = juncs.loc[juncs["strand"] == "-", "ss3"] -1
juncs["exon_end"] = "NA"
juncs.loc[juncs["strand"] == "+", "exon_end"] = juncs.loc[juncs["strand"] == "+", "ss5"] -1
juncs.loc[juncs["strand"] == "-", "exon_end"] = juncs.loc[juncs["strand"] == "-", "ss5"] +1

# collapse splice sites 
sum_cols_ss3 = ['chr',"strand", "exon_start","ss3"]
ss3_collapse= juncs.groupby("id_ss3")[sum_cols_ss3].agg(lambda x: "_".join(map(str, x.unique()))).reset_index()
sum_cols_ss5 = ['chr',"strand", "exon_end","ss5"]
ss5_collapse= juncs.groupby("id_ss5")[sum_cols_ss5].agg(lambda x: "_".join(map(str, x.unique()))).reset_index()
# convert to numeric 
ss3_collapse['exon_start'] = pd.to_numeric(ss3_collapse['exon_start'])
ss5_collapse['exon_end'] = pd.to_numeric(ss5_collapse['exon_end'])

# add bed starts and bed ends 
# for collapsed ss3
ss3_collapse["bed_start"] = "NA"
ss3_collapse["bed_end"] = "NA"
# plus strand
ss3_collapse.loc[ss3_collapse["strand"] == "+", "bed_start"] = \
    ss3_collapse.loc[ss3_collapse["strand"] == "+", "exon_start"] -1 - intron_window
ss3_collapse.loc[ss3_collapse["strand"] == "+", "bed_end"] = \
    ss3_collapse.loc[ss3_collapse["strand"] == "+", "exon_start"] + exon_window
# minus strand 
ss3_collapse.loc[ss3_collapse["strand"] == "-", "bed_start"] = \
    ss3_collapse.loc[ss3_collapse["strand"] == "-", "exon_start"] -1 - exon_window
ss3_collapse.loc[ss3_collapse["strand"] == "-", "bed_end"] = \
    ss3_collapse.loc[ss3_collapse["strand"] == "-", "exon_start"] + intron_window
# for collapsed ss5
ss5_collapse["bed_start"] = "NA"
ss5_collapse["bed_end"] = "NA"
# plus strand
ss5_collapse.loc[ss5_collapse["strand"] == "+", "bed_start"] = \
    ss5_collapse.loc[ss5_collapse["strand"] == "+", "exon_end"] -1 - exon_window
ss5_collapse.loc[ss5_collapse["strand"] == "+", "bed_end"] = \
    ss5_collapse.loc[ss5_collapse["strand"] == "+", "exon_end"] + intron_window
# minus strand 
ss5_collapse.loc[ss5_collapse["strand"] == "-", "bed_start"] = \
    ss5_collapse.loc[ss5_collapse["strand"] == "-", "exon_end"] -1 - intron_window
ss5_collapse.loc[ss5_collapse["strand"] == "-", "bed_end"] = \
    ss5_collapse.loc[ss5_collapse["strand"] == "-", "exon_end"] + exon_window

print("generating perBase coverage...")
# generate bed files 
ss3_bed = pd.DataFrame({
    'chr': ss3_collapse['chr'],
    'start': ss3_collapse['bed_start'],
    'end': ss3_collapse['bed_end'],
    'id': ss3_collapse['id_ss3'],
    'score': ".",
    'strand': ss3_collapse['strand']
})
# generate bed files 
ss5_bed = pd.DataFrame({
    'chr': ss5_collapse['chr'],
    'start': ss5_collapse['bed_start'],
    'end': ss5_collapse['bed_end'],
    'id': ss5_collapse['id_ss5'],
    'score': ".",
    'strand': ss5_collapse['strand']
})
ss3_bed_filename= outpath+"ss3.bed"
ss3_bed.to_csv(ss3_bed_filename, sep="\t", header=False, index=False)
ss5_bed_filename= outpath+"ss5.bed"
ss5_bed.to_csv(ss5_bed_filename, sep="\t", header=False, index=False)


# run coverage
# for ss3 centered
outpath_cov_ss3= outpath + "cov_ss3.txt"
cov_command_ss3 ="bedtools coverage -d -split -a "\
            + ss3_bed_filename +" -b "+ bam +" > "+outpath_cov_ss3

os.system(cov_command_ss3)
# for ss5 centered
outpath_cov_ss5= outpath + "cov_ss5.txt"
cov_command_ss5 ="bedtools coverage -d -split -a "\
            + ss5_bed_filename +" -b "+ bam +" > "+outpath_cov_ss5

os.system(cov_command_ss5)

print("predicting exon boundaries...")

# import coverage and detecting exon boundaries for ss3
cov_ss3 = pd.read_csv(outpath_cov_ss3,\
                 sep="\t",low_memory=False, header = None )  
cov_ss3.columns = ["chr","start", "end","id", "score", "strand","pos", "cov"]
cov_ss3.loc[cov_ss3["cov"] < 3, "cov"] = 0
# collapse coverage into lists
cov_collapse_ss3=(cov_ss3.iloc[:, :8]                   # first 8 columns
    .groupby(["id", "chr", "start", "end", "score","strand"], as_index=False)  # group by these columns
    .agg({"cov": list})                      # aggregate 'cov' into a list
)
def exon_predict(cov_list):
    model = "l2"  
    binseg = rpt.Binseg(model=model, jump=1).fit(np.array(cov_list))
    start_predict=binseg.predict(n_bkps=2)[0]
    end_predict=binseg.predict(n_bkps=2)[1]
    return start_predict, end_predict

#start = time.time()
cov_collapse_ss3["predict"] = [exon_predict(x) for x in cov_collapse_ss3["cov"]]
#end = time.time()
#print(end - start)
cov_collapse_ss3["start_predict"] = cov_collapse_ss3["predict"].str[0]
cov_collapse_ss3["end_predict"] = cov_collapse_ss3["predict"].str[1]

# add coordinates for predicted start and end of exons
cov_collapse_ss3["exon_predicted_start"] = cov_collapse_ss3["start_predict"] + 1 + cov_collapse_ss3['start']
cov_collapse_ss3["exon_predicted_end"] = cov_collapse_ss3["end_predict"] + cov_collapse_ss3['start']

# import coverage and detecting exon boundaries for ss5
cov_ss5 = pd.read_csv(outpath_cov_ss5,\
                 sep="\t",low_memory=False, header = None )  
cov_ss5.columns = ["chr","start", "end","id", "score", "strand","pos", "cov"]
cov_ss5.loc[cov_ss5["cov"] < 3, "cov"] = 0
# collapse coverage into lists
cov_collapse_ss5=(cov_ss5.iloc[:, :8]                   # first 8 columns
    .groupby(["id", "chr", "start", "end", "score","strand"], as_index=False)  # group by these columns
    .agg({"cov": list})                      # aggregate 'cov' into a list
)
def exon_predict(cov_list):
    model = "l2"  
    binseg = rpt.Binseg(model=model, jump=1).fit(np.array(cov_list))
    start_predict=binseg.predict(n_bkps=2)[0]
    end_predict=binseg.predict(n_bkps=2)[1]
    return start_predict, end_predict

#start = time.time()
cov_collapse_ss5["predict"] = [exon_predict(x) for x in cov_collapse_ss5["cov"]]
#end = time.time()
#print(end - start)
cov_collapse_ss5["start_predict"] = cov_collapse_ss5["predict"].str[0]
cov_collapse_ss5["end_predict"] = cov_collapse_ss5["predict"].str[1]

# add coordinates for predicted start and end of exons
cov_collapse_ss5["exon_predicted_start"] = cov_collapse_ss5["start_predict"] + cov_collapse_ss5['start'] + 1
cov_collapse_ss5["exon_predicted_end"] = cov_collapse_ss5["end_predict"] + cov_collapse_ss5['start']


# make high confidence exon list 
# make ID for ss3 centered exons 
cov_collapse_ss3['exon_id'] =  cov_collapse_ss3["chr"].astype(str) + "_" +\
                                cov_collapse_ss3["exon_predicted_start"].astype(str) + "_" +\
                                 cov_collapse_ss3["exon_predicted_end"].astype(str) + "_" +\
                                  cov_collapse_ss3["strand"].astype(str)
# get splice site coordinate for ss3
cov_collapse_ss3['centered_ss3'] = cov_collapse_ss3['id'].str.split('_').str[1].astype(int)
# make ID for ss5 centered exons 
cov_collapse_ss5['exon_id'] =  cov_collapse_ss5["chr"].astype(str) + "_" +\
                                cov_collapse_ss5["exon_predicted_start"].astype(str) + "_" +\
                                 cov_collapse_ss5["exon_predicted_end"].astype(str) + "_" +\
                                  cov_collapse_ss5["strand"].astype(str)
# get splice site coordinate for ss5
cov_collapse_ss5['centered_ss5'] = cov_collapse_ss5['id'].str.split('_').str[1].astype(int)

matched =  pd.merge(cov_collapse_ss3, cov_collapse_ss5, on='exon_id',suffixes=('_ss3', '_ss5'))

print("generating high confidence exons...")

# add expected exon start and exon end based on centered 3ss and 5ss
matched["expected_exon_start"] = "NA"
matched.loc[matched["strand_ss3"] == "+", "expected_exon_start"] = matched.loc[matched["strand_ss3"] == "+", "centered_ss3"] +1
matched.loc[matched["strand_ss3"] == "-", "expected_exon_start"] = matched.loc[matched["strand_ss3"] == "-", "centered_ss3"] -1
matched["expected_exon_end"] = "NA"
matched.loc[matched["strand_ss3"] == "+", "expected_exon_end"] = matched.loc[matched["strand_ss3"] == "+", "centered_ss5"] -1
matched.loc[matched["strand_ss3"] == "-", "expected_exon_end"] = matched.loc[matched["strand_ss3"] == "-", "centered_ss5"] +1

# add metrics for distances of expected exon boundaries based on centered splice site to predicted exon boundaries
# distance of predicted exon boundaries to expected exon boundaries based on splicesites
matched['splice_distance_ss3'] = matched['exon_predicted_start_ss3'] - matched['expected_exon_start']
matched['splice_distance_ss5'] = matched['exon_predicted_end_ss5'] - matched['expected_exon_end']
# distance between predicted exon start from  centered ss3 vs centered ss5
matched['exon_start_distance'] = matched['exon_predicted_start_ss3'] - matched['exon_predicted_start_ss5']
matched['exon_end_distance'] = matched['exon_predicted_end_ss3'] - matched['exon_predicted_end_ss5']

matched_fil = matched[(matched['splice_distance_ss3'] == 0) & (matched['splice_distance_ss5'] == 0)]

# format dataframes for output 
# formatting high confidence exons
hiconf_exons = matched_fil[['chr_ss3','centered_ss3', 'centered_ss5','strand_ss3',\
            'exon_predicted_start_ss3','exon_predicted_end_ss3'] ]

hiconf_exons.columns = ["chr","centered_ss3", "centered_ss3","strand", "predicted_exonStart", "predicted_exonEnd"]
# formating ss3 centered exons 
ss3_exons = cov_collapse_ss3[['chr','centered_ss3','strand','exon_predicted_start','exon_predicted_end']]
ss3_exons.columns = ["chr","centered_ss3", "strand", "predicted_exonStart", "predicted_exonEnd"]
# formating ss5 centered exons 
ss5_exons = cov_collapse_ss5[['chr','centered_ss5','strand','exon_predicted_start','exon_predicted_end']]
ss5_exons.columns = ["chr","centered_ss5", "strand", "predicted_exonStart", "predicted_exonEnd"]

print("writing output...") 

# generate output files 
# high confidence exons
hiconfExons_filename = outpath + "predictedExons_hiConf.txt"
hiconf_exons.to_csv(hiconfExons_filename, index=False, sep = "\t")
# exons from centered ss3
ss3_exons_filename = outpath + "predictedExons_ss3.txt"
ss3_exons.to_csv(ss3_exons_filename, index=False, sep = "\t")
# exons from centered ss5
ss5_exons_filename = outpath + "predictedExons_ss5.txt"
ss5_exons.to_csv(ss5_exons_filename, index=False, sep = "\t")

