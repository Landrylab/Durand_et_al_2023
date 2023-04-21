import pandas as pd
import numpy as np
import io
import os
import gzip
import sys



def readVcfGzip(filepath):
    fh = gzip.open(filepath,'rt')
    lines = [l for l in fh if not l.startswith('##')]
    df_vcf = pd.read_csv(io.StringIO(''.join(lines)),sep='\t')
    return(df_vcf)


def filterVcf(df_):
    """
    This filter masks haplotype genotypes with "." UNLESS:
    1) alleleic depth (AD) of the second most common variant within genotype is less than 4 
       sorted(AD, reverse=True)[1] <= 4
       AND
    2) allelic depth of the second most frequent variant to most frequent variant is less than 0.2
       sorted(AD, reverse=True)[1]/sorted(AD, reverse=True)[0] <= 0.2
       
    #3) allelic depth of reads for the variant that was called have at least one forward and one reverse read
    #   (ADF[int(GT)] > 0) & (ADR[int(GT)] > 0)
    """
    
    sample_cols = df_.columns[9:]

    DO = {}
    for sample in sample_cols:
        gts = df_[sample]
        # GT format: GT:PL:DP:ADF:ADR:AD:GQ
        OUT = []
        MC = []
        mixed_count = 0
        for genotype in gts:
            gt = genotype.split(":")
            GT = gt[0]
            DP = gt[2]
            AD = [int(i) for i in gt[5].split(',')] # allelic depth for each variant: ref, alt
            AD_nonref = sum([AD[0]] + AD[2:]) # sum of allelic depths for all variants except for ref
            sAD = sorted([int(i) for i in gt[5].split(',')], reverse=True) # allelic depth for each variant from most to least common
            sAD1 = sum(sAD[1:]) # sum of allelic depths for all variants except for the most common one
            ADF = [int(i) for i in gt[3].split(',')] # forward allelic depth for each variant
            ADR = [int(i) for i in gt[4].split(',')]  # reverse allelic depth for each variant 
            if (GT == '.'):
                out = genotype
            elif GT == "2":
                out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
            elif GT == "3":
                out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
            elif GT == "1":
                if int(DP)<4.:
                    out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
                else:
                    if sAD1>0: # sum of all alleles exept for the most common
                        #if (AD[1] <= 4) & (AD[1]/AD[0] <= 0.2) & (ADF[int(GT)] > 0) & (ADR[int(GT)] > 0):
                        #if (ADF[int(GT)] > 0) & (ADR[int(GT)] > 0):
                        if (sAD[1] <= 4) & (sAD1/sAD[0] <= 0.2) & (AD[0] <= 4) & (AD[1] > AD[0]):
                            out = genotype
                        else:
                            out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
                    else:
                        out = genotype
            elif GT == "0":
                if int(DP)<4.:
                    out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
                else:
                    if sAD1>0: # sum of all alleles exept for the most common
                        if (sAD[1] <= 4) & (sAD1/sAD[0] <= 0.2) & (AD[1] <= 4) & (AD[0] > AD[1]):
                            out = genotype
                        else:
                            out = ':'.join([".",gt[1],gt[2],gt[3],gt[4],gt[5],gt[6]])
                    else:
                        out = genotype
            #OUT[sample] = genotype,out
            OUT.append(out)

        DO[sample] = OUT

    dO = pd.DataFrame(DO)

    
    vcf_info = df_.iloc[:,0:9]
    vcf_info['POS'] = vcf_info['POS'].astype(str)
    vcf_info['QUAL'] = vcf_info['QUAL'].astype(str)
    vcf_filtered = pd.concat([vcf_info,dO], axis=1)
    return(vcf_filtered)

def writeVcf(path, df_body_):
    
    fh = gzip.open(path,'rt')
    header = ''.join([h for h in fh.readlines() if h.startswith('#')])
    filename_ID = path.replace('.vcf.gz','')
    
    wh = open(filename_ID+'_filterAD.vcf', 'w')
    wh.write(header)
    vcf_tab = df_body_.values.tolist()
    for snp in vcf_tab:
        f_snp = '\t'.join(snp)+'\n'
        wh.write(f_snp)
    wh.flush()
    wh.close()

def writeTab(path, df):
    filename_ID = path.replace('.vcf.gz','').split("/")[-1]
    df.to_csv(filename_ID+"_MixedSamples.tab", sep=",", header=True, index=False)

if __name__ == '__main__':
    
    #filename = '../temp/vcf_LL13_DP.vcf.gz'
    
    filename = sys.argv[1]
    df = readVcfGzip(filename)
    df_filtered = filterVcf(df)
    writeVcf(filename, df_filtered)
    
    
    