import numpy as np
import os
from Bio import SeqIO
import sys
import re

def sorted_alphanumeric(data):
    convert=lambda text:int(text) if text.isdigit() else text.lower()
    alphanum_key=lambda key: [convert(c) for c in re.split('([0-9]+)',key)]
    return sorted(data,key=alphanum_key)

output_dir=sys.argv[1]
this_set=sys.argv[2]
patient=sys.argv[3]
protein=sys.argv[4]

ref_header='reference_sequence, all zeros'
ref_seq='A'*Lin

indir=output_dir+'/reconstructed_haplotypes/'+this_set+'/'+patient+'/'+protein

outdir0=output_dir+'/MPL_format/'
if not os.path.exists(outdir0):
    os.makedirs(outdir0)

outdir1=outdir0+this_set+'/'
if not os.path.exists(outdir1):
    os.makedirs(outdir1)

outdir2=outdir1+patient+'/'
if not os.path.exists(outdir2):
    os.makedirs(outdir2)

outdir3=outdir2+protein+'/'
if not os.path.exists(outdir3):
    os.makedirs(outdir3)

for MC in range(1,121):
#for MC in [4,15,20,30,48,54,73,76,81,83,91,95,119]:
    in_dir=indir+'MC_run_'+str(MC)+'/'
    out_dir0=outdir+'p'+str(MC)+'/'
    if not os.path.exists(out_dir0):
        os.makedirs(out_dir0)
    out_dir=out_dir0+'synth'+'/'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    tp=range(1,11)
    #tp=range(1,2)
    for curr_tp in tp:
        with open(out_dir+'p'+str(MC)+'_'+prot+'_bsample1of1_t'+str(dT*(curr_tp)+1)+'.fasta', 'w') as f:
            f.write('>'+ref_header+'\n'+ref_seq+'\n')
            # if(os.path.isfile(in_dir+'MSA_'+str(curr_tp)+'/result_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'/'+'prefix_consensus_aln.fasta')):
            for record in SeqIO.parse(in_dir+'MSA_'+str(curr_tp)+'/result_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'/'+'prefix_consensus_aln.fasta',"fasta"):
                #print('writing aligned consensus')
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
                #print(header_split)
                if('ref' not in curr_header): # not consensus
                    if(len(header_split)>2): # actual read
                        modif_header='read'+header_split[1]+'_'+str(np.round(float(header_split[5]),4))
                    else:
                        modif_header='read0_1.0000'
                    f.write('>'+modif_header+'\n'+curr_seq+'\n')
                # else: # if there's only 1 recons sequences, prefix_consensus_aln is empty, so use the original prefix_consensus file
                #     for record in SeqIO.parse(in_dir+'MSA_'+str(curr_tp)+'/result_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'/'+'prefix_consensus.fasta',"fasta"):
                #         #print('writing consensus')
                #         curr_header=record.id
                #         curr_seq=str(record.seq).upper()
                #         header_split=curr_header.split('_')
                #         if(len(header_split)!=1): # not consensus
                #             modif_header='read'+header_split[1]+'_'+str(np.round(float(header_split[5]),4))
                #             f.write('>'+modif_header+'\n'+curr_seq+'\n')
        f.close()
        # else: # sequence in reads can't be reconstructed, so using consensus from reads in its place
        #     pysam.sort("-o",in_dir+'MSA_'+str(curr_tp)+'/'+'MSA_'+str(curr_tp)+'_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'_sorted.sam',in_dir+'MSA_'+str(curr_tp)+'/'+'MSA_'+str(curr_tp)+'_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'.sam')
        #     pysam.consensus("-o",in_dir+'MSA_'+str(curr_tp)+'/result_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'/'+'reads_consensus.fasta', in_dir+'MSA_'+str(curr_tp)+'/'+'MSA_'+str(curr_tp)+'_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'_sorted.sam')
        #     os.remove(in_dir+'MSA_'+str(curr_tp)+'/'+'MSA_'+str(curr_tp)+'_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'_sorted.sam')
        #     with open(out_dir+'p'+str(MC)+'_'+prot+'_bsample1of1_t'+str(dT*(curr_tp-1)+1)+'.fasta', 'w') as f:
        #         f.write('>'+ref_header+'\n'+ref_seq+'\n')
        #         for record in SeqIO.parse(in_dir+'MSA_'+str(curr_tp)+'/result_'+str(mean_iden)+'_'+str(read_len)+'_'+str(covg)+'/'+'reads_consensus.fasta',"fasta"):
        #             curr_seq=str(record.seq).upper()
        #             modif_header='read0_1.0000'
        #             f.write('>'+modif_header+'\n'+curr_seq+'\n')
        #     f.close()