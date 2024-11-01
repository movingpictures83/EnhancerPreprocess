#!/usr/bin/env python




######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################
def print_error(error_message=""):
    sys.stderr.write(error_message)
    sys.exit(1)

def get_epimark_profile(query_bed_file,epimark_id,bw_file,output_filename):
    intervals = pd.read_table(query_bed_file,
                              names=["chr","start","end","id"])

    cmd = " ".join(["/usr/bin/env",
                    "bigWigAverageOverBed",
                    bw_file,
                    query_bed_file,
                    "stdout"])

    ## The code regarding subprocess.Popen is based from oarevalo's answer on:
    ## http://stackoverflow.com/questions/16198546/get-exit-code-and-stderr-from-subprocess-call
    pipes = subprocess.Popen(shlex.split(cmd),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print_error(std_err)
    score = pd.read_table(StringIO(std_out),
                          names=["id","size","covered","sum","mean0",mark + "_" + sample]
    )
     
    #score = pd.read_table(
    #    StringIO(subprocess.check_output(shlex.split(cmd),
    #                                     universal_newlines=True)),
    #    names=["id","size","covered","sum","mean0",mark + "_" + sample]
    #)
        
    #intervals = intervals["id"]
    score = pd.merge(intervals,score,on="id")
    del intervals
    #score = score[["id",mark + "_" + sample]]
    score = score[mark + "_" + sample]
    score.to_csv(output_filename,sep="\t",header=False,index=False)
    del score

def merge_split_epimark_files(query_bed_file,epimark_file_prefix,
                              epimark_id_list,header,output_filename):
    outhandle = open(output_filename,'w')
    outhandle.write(header)
    num_line = 0
    inhandle_list = []
    for epimark_id in epimark_id_list:
        inhandle_list.append(open(epimark_file_prefix+str(epimark_id)+".txt",'r'))
    fhandle = open(query_bed_file,'r')
    for line in fhandle:
        num_line += 1
        line = line.rstrip()
        for i in range(len(inhandle_list)):
            line = line + "\t" + inhandle_list[i].readline().rstrip()    
        outhandle.write(line+"\n")
    for i in range(len(inhandle_list)):
        inhandle_list[i].close()
    outhandle.close()
    return(num_line)
   
######################################################################################
#Usage information
import argparse
import sys
######################################################################################
#Import modules
    
import multiprocessing
import subprocess
import shlex
import os
import pandas as pd


import PyPluMA
import PyIO
class EnhancerPreprocessPlugin:
 def input(self, inputfile):
  self.parameters = PyIO.readParameters(inputfile)
 def run(self):
  pass
 def output(self, outputfile):
  self.data_info_file = PyPluMA.prefix()+"/"+self.parameters["data_info_file"]
  self.query_region_file = PyPluMA.prefix()+"/"+self.parameters["query_region_file"]
  self.output_prefix = outputfile
  self.DMR_file = PyPluMA.prefix()+"/"+self.parameters["DMR_file"]
  self.num_process = 1
  if ("num_process" in self.parameters):
     self.num_process = int(self.parameters["num_process"])
  self.num_output_file = 1
  if ("num_output_file" in self.parameters):
         self.num_output_file = int(self.num_output_file)
  self.for_genome_wide_prediction = True
  if ("for_genome_wide_prediction" in self.parameters and self.parameters["for_genome_wide_prediction"] == "False"):
          self.for_genome_wide_prediction = False
  self.DMR_overlap_fraction = 1.0
  self.extending_DMR =  0

  #Check whether fraction (-f) is within range
  if self.DMR_overlap_fraction > 1.0 or self.DMR_overlap_fraction < 0.0:
    print("Error: Invalid value for -f: " + str(self.DMR_overlap_fraction))
    print("       Should be between 0.0 and 1.0")

  try:
    from StringIO import StringIO
  except:
    from io import StringIO
    
######################################################################################
#Read sample information
  try:
    data_info = pd.read_table(self.data_info_file)
  except:
    print_error("Error! Cannot open file \"" + self.data_info_file + "\"!\n")

#Check if there are duplicated sample and mark combination
  if sum(data_info.duplicated(["sample","mark"])) > 0:
    print(data_info.ix[data_info.duplicated(["sample","mark"])])
    print_error("Error! Found duplicated sample and mark combination!\n")

#Check if every bigwig file exists
  for (sample,mark,bw_file) in data_info.itertuples(index=False):
    if not (os.path.isfile(PyPluMA.prefix()+"/"+bw_file)):
        print_error("Error! bigwig file \"" + bw_file + "\" not exist!\n")

######################################################################################
  if self.DMR_file is not None:
    if not self.for_genome_wide_prediction:
        cmd = " ".join(["/usr/bin/env",
                        "bedtools",
                        "intersect",
                        "-a",self.DMR_file,
                        "-b",self.query_region_file,
                        "-wa -wb",
                        "-f",str(self.DMR_overlap_fraction)])

        #Get DMRs that overlap regions
        ## The code regarding subprocess.Popen is based from oarevalo's answer on:
        ## http://stackoverflow.com/questions/16198546/get-exit-code-and-stderr-from-subprocess-call
        pipes = subprocess.Popen(shlex.split(cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        std_out, std_err = pipes.communicate()
        if pipes.returncode != 0:
            print_error(std_err)

        DMR = pd.read_table(StringIO(std_out),
                            names = ["dmr_chr","dmr_start","dmr_end","dmr_id",
                                     "reg_chr","reg_start","reg_end","reg_id"])

        DMR["combined_id"] = DMR["dmr_id"]+":"+DMR["reg_id"]
        DMR = DMR[["dmr_chr","dmr_start","dmr_end","combined_id"]]
    else:
        DMR = pd.read_table(self.DMR_file,
                            names = ["dmr_chr","dmr_start","dmr_end","dmr_id"])

    DMR["dmr_start"] = DMR["dmr_start"] - self.extending_DMR
    DMR.loc[DMR["dmr_start"] < 0,"dmr_start"] = 0
    DMR["dmr_end"] = DMR["dmr_end"] + self.extending_DMR
    DMR.to_csv(self.output_prefix+".DMR.bed",sep="\t",header=False,index=False)
    del DMR
######################################################################################
#
#print(min(self.num_process,data_info.shape[0]))
  pool = multiprocessing.Pool(min(self.num_process,data_info.shape[0]))

  for (sample,mark,bw_file) in data_info.itertuples(index=False):
    epimark_id = mark + "_" + sample
    pool.apply_async(
        get_epimark_profile,
        (self.query_region_file,
         epimark_id,
         PyPluMA.prefix()+"/"+bw_file,
         self.output_prefix+".region."+str(epimark_id)+".txt"),
    )
    if self.DMR_file is not None:
        pool.apply_async(
            get_epimark_profile,
            (self.output_prefix+".DMR.bed",
             epimark_id,
             PyPluMA.prefix()+"/"+bw_file,
             self.output_prefix+".DMR."+str(epimark_id)+".txt"),
        )


  pool.close()
  pool.join()
  ######################################################################################
  # Merging
  ## get header
  epimark_id_list = [mark+"_"+sample for (sample,mark,bw_file) in data_info.itertuples(index=False)]
  header = "\t".join(["chr","start","end","id"]+epimark_id_list) + "\n"

## begin merging
  pool = multiprocessing.Pool(min(self.num_process,2))

  num_region = pool.apply_async(
    merge_split_epimark_files,
    (self.query_region_file,
     self.output_prefix+".region.",
     epimark_id_list,
     header,
     self.output_prefix+".region_with_epimark.tsv"),
  )

  if self.DMR_file is not None:
    num_DMR = pool.apply_async(
        merge_split_epimark_files,
        (self.output_prefix+".DMR.bed",
         self.output_prefix+".DMR.",
         epimark_id_list,
         header,
         self.output_prefix+".DMR_with_epimark.tsv"),
    )

  pool.close()
  pool.join()
  num_region = num_region.get(timeout=1)

  #if self.DMR_file is not None:
  #  num_DMR = num_DMR.get(timeout=1)
  #  subprocess.check_call(["rm"]+
  #                        [self.output_prefix+".DMR.bed"]+
  #                        [self.output_prefix+".DMR."+str(epimark_id)+".txt" for epimark_id in epimark_id_list]+
  #                        [self.output_prefix+".region."+str(epimark_id)+".txt" for epimark_id in epimark_id_list])
  #else:
  #  subprocess.check_call(["rm"]+
  #                        [self.output_prefix+".region."+str(epimark_id)+".txt" for epimark_id in epimark_id_list])

######################################################################################
# Split files
  if self.num_output_file > 1:
    out_line_cutoff = int(float(num_region)/float(self.num_output_file))+1
    out_line_count = 0
    out_file_id = 0
    fhandle = open(self.output_prefix+".region_with_epimark.tsv",'r')
    fhandle.readline()
    outhandle = open(self.output_prefix+".region_with_epimark."+str(out_file_id)+".tsv",'w')
    outhandle.write(header)
    regid2fileid = {}
    for line in fhandle:
        if out_line_count >= out_line_cutoff:
            outhandle.close()
            out_line_count = 0
            out_file_id += 1
            outhandle = open(self.output_prefix+".region_with_epimark."+str(out_file_id)+".tsv",'w')
            outhandle.write(header)
        outhandle.write(line)
        out_line_count += 1
        regid2fileid[line.split("\t")[3]] = out_file_id
        
    outhandle.close()
    fhandle.close()        

    ## DMR
    if self.DMR_file is not None:
        if not self.for_genome_wide_prediction:
            fhandle = open(self.output_prefix+".DMR_with_epimark.tsv",'r')
            fhandle.readline() ## Skip header
            out_line_count = 0
            outhandles = []
            for out_file_id in range(self.num_output_file):
                outhandles.append(open(self.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w'))
                outhandles[-1].write(header)
        
            for line in fhandle:
                region_id = line.split("\t")[3].split(":")[1]
                outhandles[regid2fileid[region_id]].write(line)
            fhandle.close()
    
            for i in range(len(outhandles)):
                outhandles[i].close()

        else:
            ## For genome-wide enhancer prediction
            ## DMR (without link to query regions)
            out_line_cutoff = int(float(num_DMR)/float(self.num_output_file))+1
            out_line_count = 0
            out_file_id = 0
            fhandle = open(self.output_prefix+".DMR_with_epimark.tsv",'r')
            fhandle.readline()
            outhandle = open(self.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w')
            outhandle.write(header)
            regid2fileid = {}
            for line in fhandle:
                if out_line_count >= out_line_cutoff:
                    outhandle.close()
                    out_line_count = 0
                    out_file_id += 1
                    outhandle = open(self.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w')
                    outhandle.write(header)
                outhandle.write(line)
                out_line_count += 1
            
            outhandle.close()
            fhandle.close()
