import os
import sys
import re

def sorted_alphanumeric(data):
    convert=lambda text:int(text) if text.isdigit() else text.lower()
    alphanum_key=lambda key: [convert(c) for c in re.split('([0-9]+)',key)]
    return sorted(data,key=alphanum_key)

data_dir=sys.argv[1]
output_dir=sys.argv[2]
bash_dir=sys.argv[3]
logs_dir=sys.argv[4]
rvhaplo_dir=sys.argv[5]
this_set=sys.argv[6]
patient=sys.argv[7]
protein=sys.argv[8]
numcores=int(sys.argv[9])

mainoutdir_main=bash_dir+'reconstruction_scripts'+'/'
if(not os.path.isdir(mainoutdir_main)):    
    os.mkdir(mainoutdir_main)
mainoutdir0=mainoutdir_main+this_set+'/'
if(not os.path.isdir(mainoutdir0)):
    os.mkdir(mainoutdir0)            
mainoutdir1=mainoutdir0+patient+'/'
if(not os.path.isdir(mainoutdir1)):    
    os.mkdir(mainoutdir1)
mainoutdir=mainoutdir1+protein+'/'
if(not os.path.isdir(mainoutdir)):
    os.mkdir(mainoutdir)

logsoutdir0=logs_dir+this_set+'/'
if(not os.path.isdir(logsoutdir0)):    
    os.mkdir(logsoutdir0)
logsoutdir1=logsoutdir0+patient+'/'
if(not os.path.isdir(logsoutdir1)):    
    os.mkdir(logsoutdir1)              
logsoutdir=logsoutdir1+protein+'/'
if(not os.path.isdir(logsoutdir)):    
    os.mkdir(logsoutdir)     

reconsdir_main=output_dir+'reconstructed_haplotypes'+'/'
if(not os.path.isdir(reconsdir_main)):    
    os.mkdir(reconsdir_main)
reconsdir0=reconsdir_main+this_set+'/'
if(not os.path.isdir(reconsdir0)):
    os.mkdir(reconsdir0)            
reconsdir1=reconsdir0+patient+'/'
if(not os.path.isdir(reconsdir1)):    
    os.mkdir(reconsdir1)
mainreconsdir=reconsdir1+protein+'/'
if(not os.path.isdir(mainreconsdir)):
    os.mkdir(mainreconsdir)         

mainoutdir_call0=bash_dir+'reconstruction_call'+'/'
if(not os.path.isdir(mainoutdir_call0)):    
    os.mkdir(mainoutdir_call0)
mainoutdir_call=mainoutdir_call0+this_set
if(not os.path.isdir(mainoutdir_call)):    
    os.mkdir(mainoutdir_call)
    
ref_dir=data_dir+'reads'+'/'+this_set+'/'+'ref'+'/'+patient+'_'+'ref'+'/'+protein+'/'    
ref_file=ref_dir+patient+'_'+protein+'_'+'ref'+'.fa'
maindir=output_dir+'BAM'+'/'+this_set+'/'+patient+'/'+protein
mainfile_name=mainoutdir_call+'/'+this_set+'_'+patient+'_'+protein+'_'+'reconstruction_call.sh'
filemain = open(os.path.join(mainfile_name),'w',newline='\n')
filemain.write('#!/bin/sh')
filemain.write('\n\n')
filemain.write('chmod -R 700 '+mainoutdir)
filemain.write('\n\n')
fname_count=0
pat_tps=sorted_alphanumeric(os.listdir(maindir))
for dirname in pat_tps:
    source_dir=maindir+'/'+dirname+'/'
    source_file=source_dir+dirname+'_sorted.sam'
    dest_file=mainoutdir+str(fname_count+1)+'.sh'
    dest_file_log=logsoutdir+dirname    
    fileID = open(os.path.join(dest_file),'w',newline='\n')
    fileID.write('#!/bin/sh')
    fileID.write('\n\n')   
    fileID.write('cd "'+rvhaplo_dir+'"')
    fileID.write('\n')
    fileID.write('{ time ./rvhaplo.sh -i "'+source_file+'" -r "'+ref_file+'" -o "'+mainreconsdir+dirname+'" -p '+dirname+' -t '+str(numcores)+' -wr 0.85 -wc 0.85 -ss 1 -l 0 >/dev/null 2>&1; } 2>> "'+dest_file_log+'_recons_time.txt"')
    fileID.write('\n')
    fileID.write('exit 0')
    fileID.close()
    filemain.write(os.path.join(dest_file)+' ')
    filemain.write('\n')
    fname_count=fname_count+1
filemain.write('\n')
filemain.write('exit 0')
filemain.close()