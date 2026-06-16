import numpy as np
import os
from Bio import SeqIO
import re
import sys

def sorted_alphanumeric(data):
    convert=lambda text:int(text) if text.isdigit() else text.lower()
    alphanum_key=lambda key: [convert(c) for c in re.split('([0-9]+)',key)]
    return sorted(data,key=alphanum_key)

output_dir=sys.argv[1]
this_set=sys.argv[2]
patient=sys.argv[3]
protein=sys.argv[4]

ref_header='reference_sequence, all zeros'

indir=output_dir+'aligned_haplotypes/'+this_set+'/'+patient+'/'+protein

outdir0=output_dir+'MPL_format/'
if not os.path.exists(outdir0):
    os.makedirs(outdir0)

outdir1=outdir0+this_set+'/'
if not os.path.exists(outdir1):
    os.makedirs(outdir1)

outdir2=outdir1+patient+'/'
if not os.path.exists(outdir2):
    os.makedirs(outdir2)

outdir=outdir2+protein+'/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
pat_tps=sorted_alphanumeric(os.listdir(indir))
for dirname in pat_tps:
    source_dir=indir+'/'+dirname+'/'
    source_file=source_dir+dirname+'_consensus_aln.fasta'
    fname='_'.join(dirname.split('_')[0:-1])+'_bsample1of1_'+dirname.split('_')[-1]
    dest_file=outdir+fname+'.fasta'
    
    with open(dest_file, 'w') as f:
        for record in SeqIO.parse(source_file,"fasta"):
            curr_header=record.id
            curr_seq=str(record.seq).upper()            
            curr_seq=curr_seq.replace('N','-')
            curr_seq=curr_seq.replace('R','-')
            curr_seq=curr_seq.replace('Y','-')
            curr_seq=curr_seq.replace('K','-')
            curr_seq=curr_seq.replace('M','-')
            curr_seq=curr_seq.replace('S','-')                
            curr_seq=curr_seq.replace('W','-')
            curr_seq=curr_seq.replace('B','-')
            curr_seq=curr_seq.replace('D','-')
            curr_seq=curr_seq.replace('H','-')
            curr_seq=curr_seq.replace('V','-')            
            header_split=curr_header.split('_')
            if('ref' in curr_header): # not consensus
                Lin=len(curr_seq)
                ref_seq=curr_seq
                f.write('>'+ref_header+'\n'+'A'*Lin+'\n')            
            else:
                if(len(header_split)>2): # actual read
                        modif_header='read'+header_split[1]+'_'+str(np.round(float(header_split[5]),4))
                else:
                    modif_header='read0_1.0000'
                curr_seq=''.join(cf if cg=='-' else cg for cg,cf in zip(curr_seq,ref_seq))
                f.write('>'+modif_header+'\n'+curr_seq+'\n')
    f.close()