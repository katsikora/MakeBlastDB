import os
import sys
import time
import argparse
import re
import logging
import subprocess
import fnmatch

sys.path.insert(0, "/data/boehm/group/pipelines/ruffus")
from ruffus import *
from ruffus.proxy_logger import *
from ruffus.combinatorics import *
from ruffus.drmaa_wrapper import run_job, run_job_using_drmaa, error_drmaa_job

from PIL import Image
import string


parser = argparse.ArgumentParser(prog="BlastDBpipe", version="0.0.1", description="makes a nucleotide blast DB from PE reads", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ri", "--readIn", dest="readdir", action="store", default=False, help="input read folder")
parser.add_argument("-wd", "--wdir", dest="wdir", action="store", default=False, help="output folder")
parser.add_argument("-bs", "--batchSize", dest="bsize", action="store",metavar="INT",type=int,default=10, help="number of samples to process in parallel")
parser.add_argument("-nt", "--numThr", dest="nthreads", action="store",metavar="INT",type=int,default=8, help="number of threads to use per sample")
parser.add_argument("--touchOnly", dest="touchOnly", action="store_true", help="only touch files")
parser.add_argument("--target_tasks", dest="target_tasks", action="store",default=[], help="target tasks")
parser.add_argument("--forcedtorun_tasks", dest="forcedtorun_tasks", action="store",default=[], help="forced to run tasks")
args = parser.parse_args()
#args = parser.parse_args(['--readIn','/data/boehm/sequencing_data/180109_SN7001180_0346_AH527VBCX2/Project_239_Holland_Boehm','--wdir','/data/processing/sikora/holland/A239','--touchOnly'])
#args


#setup central working directory
wdir=args.wdir
if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)


#setup logging
logger = logging.getLogger(__name__)
fhandler = logging.FileHandler(filename=os.path.join(wdir,'pipeline.log'), mode='a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.debug(subprocess.check_output('echo $DRMAA_LIBRARY_PATH',shell=True))



import drmaa

#initiate 1 drmaa session for the whole pipeline
mySession=drmaa.Session()
mySession.initialize()


#identify pipeline input files
readdir=args.readdir
#os.chdir(readdir)

libset2 = []
# Walk through directory
for dName, sdName, fList in os.walk(readdir):
    for fileName in fList:
        if fnmatch.fnmatch(fileName, "*fastq.gz"): # Match search string
            libset2.append(os.path.join(dName, fileName))

libset2_R1=filter(lambda x:'_R1.fastq.gz' in x, libset2)
libset2_R1.sort()
libset2_R2=filter(lambda x:'_R2.fastq.gz' in x, libset2)
libset2_R2.sort()
read_root=[ re.sub('_R1.fastq.gz','',R1f) for R1f in libset2_R1 ]
INfiles=list(zip(libset2_R1,libset2_R2))
    
logger.debug(INfiles)    


##################PATHS TO EXECUTABLES###############################################################
FQCpath='/package/FastQC-0.11.3'
cutpath='/package/cutadapt-1.9.1/bin'
prinpath='/data/boehm/sikora/tools/prinseq-lite-0.20.4' #/data/boehm/group/pipelines/tools/COMPLETE_PATH
flashpath='/data/boehm/sikora/tools/FLASH-1.2.11'
blastpath='/data/boehm/sikora/tools/ncbi-blast-2.2.30+/bin'


########adapter-trim and BQ-trim sequencing reads########
cutout=os.path.join(wdir,'reads_cut')
fqcout=os.path.join(wdir,'fastqc_cut')
@mkdir(cutout,os.path.join(cutout,'logs'))
@transform(INfiles,suffix('_R1.fastq.gz'),['_R1.fastq.gz','_R2.fastq.gz'],output_dir=cutout)
def cut_reads(infiles, outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(cutpath,'cutadapt') + ' -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5  -o ' + oo1 + ' -p ' + oo2 + ' ' + ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(cutout,"logs","%s.cut_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.cut_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = bshcmd,
                                      job_name          = 'cut_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("Cut_reads error: %s" % err)
            raise
        else:
           logger.info('Adapter trimming complete')


@transform(cut_reads,suffix('_R1.fastq.gz'),['_prin_1.fastq.gz','_prin_2.fastq.gz'],output_dir=cutout)
def trim_BQ(infiles,outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    uzcmd1='zcat -v '+ ii1 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1)))
    uzcmd2='zcat -v '+ ii2 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2)))
    bshcmd='perl '+ os.path.join(prinpath,'prinseq-lite.pl') + ' -fastq ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' -fastq2 ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' -out_good ' + os.path.join(cutout,re.sub('_1.fastq.gz','',os.path.basename(oo1))) +' -trim_qual_right 20 -trim_qual_type min -trim_qual_window 6 -trim_qual_step 3 -min_len 50 -ns_max_p 10 -min_qual_mean 26 -out_bad null'
    zcmd1='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' > ' + oo1
    zcmd2='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2))) + ' > ' + oo2
    clcmd='rm -v '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2)))
    cmd_all=';'.join([uzcmd1,uzcmd2,bshcmd,zcmd1,zcmd2,clcmd])
    logger.info(cmd_all)           
    with open(os.path.join(cutout,"logs","%s.BQtrim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.BQtrim_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = cmd_all,
                                      job_name          = 'BQtrim_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("BQtrim_reads error: %s" % err)
            raise
        else:
           logger.info('Base quality trimming complete')


@mkdir(fqcout,os.path.join(fqcout,'logs'))            
@transform(trim_BQ,suffix('_1.fastq.gz'),output=['_1.zip','_1.html','_2.zip','_2.html'],output_dir=fqcout)
def postTrim_fqc(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(FQCpath,'fastqc ')+' --outdir ' + fqcout + ' -t 8 '+ ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(fqcout,"logs","%s.post_fqc.out" % read_root),'w+') as stdoutF, open(os.path.join(fqcout,"logs","%s.post_fqc.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'post_fqc',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Post_trim_fastqc error: %s" % err)
            raise
        else:
            logger.info('Post trim fastqc complete')  


#######merge forward and reverse mates########
@transform(trim_BQ,suffix('_1.fastq.gz'),output=['_flash.extendedFrags.fastq.gz','_flash.notCombined_1.fastq.gz','_flash.notCombined_2.fastq.gz'],output_dir=cutout)
def merge_mates(input_files,output_file):
    ii1 = input_files[0]
    ii2 = input_files[1]
    oo = re.sub('.extendedFrags.fastq.gz','',os.path.basename(output_file[0]))
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(flashpath,'flash')+ ' -z -M 300 -t 8 -o '+ oo + ' -d ' + cutout + ' ' + ii1 + ' ' + ii2 
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.flash.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.flash.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'flash',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Flash error: %s" % err)
            raise
        else:
            logger.info('Merging mates complete')


####modify read names to replace spaces " " with underscores "_"
@transform(merge_mates,suffix('_prin_flash.extendedFrags.fastq.gz'),output=['_prin_flash.extendedFrags.sed.fastq.gz','_prin_flash.notCombined_1.sed.fastq.gz','_prin_flash.notCombined_2.sed.fastq.gz'],output_dir=cutout)#
def mod_Rnames(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    ii3 = input_files[2]
    oo1 = output_files[0]
    oo2 = output_files[1]
    oo3 = output_files[2]
    read_root=re.sub('_prin_flash.extendedFrags.fastq.gz','',os.path.basename(ii1))
    cmd1='zcat ' + ii1 + ' | sed \'s/\ /_/g\' - | gzip -c  > ' + oo1
    cmd2='zcat ' + ii2 + ' | sed \'s/\ /_/g\' - | gzip -c  > ' + oo2
    cmd3='zcat ' + ii3 + ' | sed \'s/\ /_/g\' - | gzip -c  > ' + oo3
    cmd_all=[cmd1,cmd2,cmd3]
    bshcmd=' ; '.join(cmd_all)
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.sed.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.sed.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'sed',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Sed error: %s" % err)
            raise
        else:
            logger.info('Renaming reads complete')


###make blast database
DBout=os.path.join(wdir,'blastDB')
if not os.path.exists(DBout):
    os.makedirs(DBout)
os.chdir(DBout)
@follows(mkdir(DBout),mkdir(os.path.join(DBout,'logs')))
@transform(mod_Rnames,suffix("_prin_flash.extendedFrags.sed.fastq.gz"),output=".flash.db.nal",output_dir=DBout)
def makeDB(input_files,output_file):
    ii1 = input_files[0]
    ii2 = input_files[1]
    ii3 = input_files[2]
    oo = output_file
    read_root=re.sub('_prin_flash.extendedFrags.sed.fastq.gz','',os.path.basename(ii1))
    oox=os.path.join(os.path.dirname(oo),(read_root + '.readDB'),read_root + '.flash.db')
    bshcmd='zcat -v ' + ii1 + ' ' + ii2 + ' ' + ii3 + ' | awk \'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}\' - | ' + os.path.join(blastpath,'makeblastdb ') + ' -in - -parse_seqids -dbtype nucl -out ' + oox + ' -title ' + read_root + '; ln -fs ' + oox + '.nal ' + oo
    logger.info(bshcmd)
    with open(os.path.join(DBout,"logs","%s.makeDB.out" % read_root),'w+') as stdoutF, open(os.path.join(DBout,"logs","%s.makeDB.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'makeDB',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("MakeDB error: %s" % err)
            raise
        else:
            logger.info('Making database complete')



#####main

if __name__ == '__main__':
    with open(os.path.join(wdir,"pipelineGraph.png"),'w') as pipeGraph:
        pipeline_printout_graph(stream=pipeGraph,output_format='png',pipeline_name='WGBS',target_tasks=args.target_tasks)
    with open (os.path.join(wdir,"pipelinePrint.txt"),'w') as pipePrint:
        pipeline_printout(verbose_abbreviated_path=0,output_stream=pipePrint,target_tasks=args.target_tasks)    

    pipeline_run(touch_files_only=args.touchOnly,multiprocess=args.bsize,target_tasks=args.target_tasks,forcedtorun_tasks=args.forcedtorun_tasks,logger=logger)
    mySession.exit()  




