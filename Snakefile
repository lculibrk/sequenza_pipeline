import re


import glob
import os
import csv
import pandas
import yaml
import re

genfile = "/projects/lculibrk_prj/CNV/Ploidetect/hg19.genome.nochr"
#configfile: "runcoverage_config_allPOGs.yaml"
windowlevels=[
        1
        ]


POGdatadir=glob.glob("/projects/POG/POG_data/POG[0-9]*[0-9]/")
pedsPOGs=glob.glob("/projects/peds_POG/POG[0-9]*[0-9]/")

POGdatadir.extend(pedsPOGs)

POGdatadir = POGdatadir

#print(pedsPOGs)

#POGflatfiles=[glob.glob(os.path.join(dir, "flatfile/*")) for dir in POGdatadir]
#POGflatfiles=[flatfile[0] for flatfile in POGflatfiles if flatfile]

Libraries = {}
Comparisons = {}
lohpaths = {}
for data_dir in POGdatadir:
    pogcase=glob.glob(os.path.join(data_dir, "flatfile/*"))
#    print(pogcase)
    if len(pogcase) > 0:
        pogcase = pogcase[0]
    else:
        continue
#    print(pogcase)
    ## Get POG id
    case = os.path.basename(re.sub(".tab", "", pogcase))
    ## Open flatfile with pandas
    with open(pogcase) as handle:
         file = pandas.read_csv(handle, sep = "\t")
    ## Grab the merged_bam entries
    file = file[file["merged_bam"].notna()]
    ## Grab lines that correspond to "WGS" and "Diseased" samples
    tumour_libs = file["library_name"][(file["diseased_status"] == "Diseased") & (file["protocol"].str.contains("WGS"))].values.tolist()
    ## Initialize empty list
    loh_tumourlibs = []
    ## Iterate over all tumour_libs
    for lib in tumour_libs:
    ## Get directory for LOH data
        #lohdir = glob.glob(os.path.join("/projects/POG/POG_data", case, "wgs", "*t_" + lib + "*n*/loh"))
#        print(lohdir)
        lohdir = glob.glob(os.path.join(data_dir, "wgs", "*t_" + lib + "*n*/loh"))
## If LOH data doesn't exist, skip this case
        if len(lohdir) == 0:
#            print("skipped! (no LOH data)")
            continue
        ## Append library to the list
        loh_tumourlibs.append(lib)
    tumour_libs = loh_tumourlibs
#    print(tumour_libs)
    ## Get the name of the normal library
    normal_lib = file["diseased_status"] == "Normal"
    normal_lib = file["library_name"][normal_lib].values.tolist()
    ## Skip case if no normal
    if len(normal_lib) == 0:
#       print("skipped! (no normal)")
       continue
    ## Make a copy of tumour_libs
    all_libs = tumour_libs.copy()
    ## Add normal libs to the list
    all_libs.extend(normal_lib)
#    print(all_libs)
    ## Make list of all tumour normal pairs
    case_comparisons = expand("{tumour}_{normal}", tumour = tumour_libs, normal = normal_lib)
    for comp in case_comparisons:
        ## Get the libs from the pair
        comparison = comp.split("_")
        ## Get the directory for LOH data
        lohdir = glob.glob(os.path.join(data_dir, "wgs", "*t_" + comparison[0] + "*n_" + comparison[1] + "*/loh"))
#        print(data_dir)
        if len(lohdir) == 0:
	## Double check that this tumour:normal pair was actually analyzed
#            print("skipped! (invalid comparison)")
            continue
        ## Add to list
        Comparisons[comp] = case
        ## Add LOH dir to dictionary of loh dirs
        lohpaths[comparison[0]] = lohdir[0]
    ## Create dictionary of bam files
    for lib in all_libs:
        Libraries[lib] = {"POG":case, "state":file["diseased_status"][(file["library_name"] == lib)].values[0], "dir":file["merged_bam"][(file["library_name"] == lib)].values[0], "library":lib}

#print(lohpaths)
#coveragelist=[]
#for library in config["POG_List"]:
#    coveragelist.extend(expand("coverage/{POG}/binned/{windowlevels}/{library}/{library}.bed", library=library, POG=config["Libraries"][library]["POG"], windowlevels=windowlevels))

#combbedlist=[]
#lohlist=[]
#lohinput=[]
#models=[]
#reports=[]
#cnas=[]
#comparisons=[]
#for comparison in Comparisons.keys():
#    if comparison == "P00169_P00129":
#        continue
#    pog=Comparisons[comparison]
#    #comparisons.extend(expand("coverage/{POG}/binned/1/{comp}_processed.bed", POG=pog, comp=comparison))    
#    #combbedlist.extend(expand("coverage/{POG}/binned/50000/{comp}_processed.bed", POG=pog, comp = comparison))
#    models.extend(expand("Ploidetect_testbuild_outs/{POG}/{comp}/{POG}-{comp}_models.txt", POG=pog, comp = comparison))
#    reports.extend(expand("Ploidetect_testbuild_outs/{POG}/{comp}/{POG}-{comp}_report.pdf", POG=pog, comp = comparison))
#    cnas.extend(expand("Ploidetect_testbuild_outs/{POG}/{comp}/{POG}-{comp}_cna_genes.txt", POG=pog, comp = comparison))

seqzlist=[]
for comparison in Comparisons.keys():
    if comparison == "P00169_P00129":
        continue
    pog=Comparisons[comparison]
#    seqzlist.extend(expand("data/seqz/{POG}/{comp}.binned.seqz.gz", POG = pog, comp = comparison))
    seqzlist.extend(expand("results/{POG}/{comp}", POG = pog, comp = comparison))


print(Comparisons)
print(Libraries)

#germline=[]
#somatic=[]
#for library in config["Libraries"]:
#    pog = config["Libraries"][library]["POG"]
#    lib = config["Libraries"][library]["library"]
#    if(config["Libraries"][library]["state"] == "Normal"):
#        germline.extend(expand("coverage/{POG}/binned/{windowlevels}/{lib}.bed", POG=pog, lib=lib, windowlevels = windowlevels))
#    if(config["Libraries"][library]["state"] == "Diseased"):
#        somatic.extend(expand("coverage/{POG}/binned/{windowlevels}/{lib}.bed", POG=pog, lib=lib, windowlevels = windowlevels))




dummy = ".gitignore"

chromosomes=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]

rule all:
    input:
        seqzlist
        #expand("data/mpileup/{lib}.mpileup.gz", lib = Libraries.keys())

rule seqz:
    input:
        lambda wildcards: Libraries[wildcards.normlib]["dir"],
        lambda wildcards: Libraries[wildcards.tumrlib]["dir"]
    output:
        "data/seqz/{pogid}/{tumrlib}_{normlib}.binned.seqz.gz"
    params:
        gc="resources/hg19gc50.txt.gz",
        fasta="resources/hg19.fa"
    conda:
        "conda_config/conda.yaml"
    shell:
        "sequenza-utils seqz_binning -s <(sequenza-utils bam2seqz -F {params.fasta} -gc {params.gc} -n {input[0]}/*dupsFlagged.bam -t {input[1]}/*dupsFlagged.bam) -w 250 | gzip > {output} "

#rule binseqz:
#    input:
#        "data/seqz/{pogid}/{tumrlib}_{normlib}.seqz.gz"
#    output:
#        "data/seqz/{pogid}/{tumrlib}_{normlib}.binned.seqz.gz"
#    shell:
#        "/projects/lculibrk_prj/CNV/sequenza/venv/bin/sequenza-utils seqz_binning -s {input} -w 250 | gzip > {output}"

rule sequenza:
    input:
        "data/seqz/{pogid}/{tumrlib}_{normlib}.binned.seqz.gz"
    output:
        directory("results/{pogid}/{tumrlib}_{normlib}")
    conda:
      	"conda_config/conda.yaml"
    benchmark:
        "benchmark/{pogid}_{tumrlib}_{normlib}.txt"
    shell:
        "lsb_release -a && Rscript scripts/run_sequenza.R -i {input} -p {wildcards.tumrlib}_{wildcards.normlib} -o {output}"
rule dummy:
    input:
    output:
        ".gitignore"
