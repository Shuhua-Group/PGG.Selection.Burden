import sys
import gzip
from collections import defaultdict
import re

def BurdenTest(vcf,damage,kegg,prefix):
    with open(damage, 'r') as f:
        snplist = defaultdict(list)
        score = {}
        header=['chr','pos','Ancestral_allele','ref','alt','clinvar_clnsig','SIFT_pred','Polyphen2_HDIV_pred','Polyphen2_HVAR_pred','CADD_phred','GERP_RS','Aloft_pred','Aloft_Confidence']
        for line in f:
            line = line.strip().split('\t')
            snp = ":".join([line[i] for i in [0,1,3,4]])
            snplist[snp] = ['sum']
            score[snp] = float(line[header.index('CADD_phred')])
            if "Likely_pathogenic" in line[header.index('clinvar_clnsig')] or "Pathogenic" in line[header.index('clinvar_clnsig')]:
                snplist[snp].append('clinvar')
            if "D" in line[header.index('SIFT_pred')]:
                snplist[snp].append('SIFT')
            if "D" in line[header.index('Polyphen2_HDIV_pred')] or "P" in line[header.index('Polyphen2_HDIV_pred')] or "D" in line[header.index('Polyphen2_HVAR_pred')] or "P" in line[header.index('Polyphen2_HVAR_pred')]:
                snplist[snp].append('Polyphen2')
            if float(line[header.index('CADD_phred')]) >= 10:
                snplist[snp].append('CADD')
                if float(line[header.index('CADD_phred')]) >= 20: snplist[snp].append('CADD20')
                else : snplist[snp].append('CADD10')
            if line[header.index('GERP_RS')] != "." and float(line[header.index('GERP_RS')]) >= 2:
                snplist[snp].append('GERP')
                if float(line[header.index('GERP_RS')]) >= 4: snplist[snp].append('GERP4')
                else : snplist[snp].append('GERP2')
            for i in range(0,len(line[header.index('Aloft_pred')].split(';'))):
                if line[header.index('Aloft_pred')].split(';')[i]=="Recessive" and line[header.index('Aloft_Confidence')].split(';')[i]=="High":
                    if 'Aloft_Recessive' not in snplist[snp]: snplist[snp].append('Aloft_Recessive')
                if line[header.index('Aloft_pred')].split(';')[i]=="Dominant" and line[header.index('Aloft_Confidence')].split(';')[i]=="High":
                    if 'Aloft_Dominant' not in snplist[snp]: snplist[snp].append('Aloft_Dominant')
    groups = {group: defaultdict(list) for group in ['sum','clinvar','SIFT','Polyphen2','CADD','CADD10','CADD20','GERP','GERP2','GERP4','Aloft_Recessive','Aloft_Dominant']} 
    CADD_weighted_groups = {group: defaultdict(list) for group in ['sum','clinvar','SIFT','Polyphen2','CADD','CADD10','CADD20','GERP','GERP2','GERP4','Aloft_Recessive','Aloft_Dominant']} 


    with open(kegg, 'r') as f:
        # kegg = defaultdict(dict)
        # genes = []
        # pathways = []
        genes = defaultdict(list)
        pathways = defaultdict(list)
        for line in f:
            line = line.strip().split('\t')
            # if line[4] not in kegg: kegg[line[4]]={}
            # kegg[line[4]].update({line[0]:[line[1],int(line[2]),int(line[3])]})
            # if line[0] not in genes: genes.append(line[0])
            # if line[4] not in pathways: pathways.append(line[4])
            if line[0] not in genes: genes[line[0]]=[line[1],int(line[2]),int(line[3])]
            if line[4] not in pathways: pathways[line[4]]=[] 
            pathways[line[4]].append(line[0])
    kegg_gene = {gene: defaultdict(list) for gene in genes}        
    kegg_pathway = {pathway: defaultdict(list) for pathway in pathways}
    total_gene = {gene: 0 for gene in genes}        
    total_pathway = {pathway: 0 for pathway in pathways}  


    with gzip.open(vcf, 'rb') as f:
        total=0
        lof=0
        for line in f:
            line = line.decode('utf-8')
            if line.startswith('##'): continue
            line = line.strip().split('\t')
            if line[0].startswith('#'):
                header = line; samples = header[9:]
                for group in groups: groups[group]=[0 for i in range(0,2*len(samples))]
                for group in CADD_weighted_groups: CADD_weighted_groups[group]=[0 for i in range(0,2*len(samples))]
                for gene in kegg_gene: kegg_gene[gene]=[0 for i in range(0,2*len(samples))] 
                for pathway in kegg_pathway: kegg_pathway[pathway]=[0 for i in range(0,2*len(samples))]
                continue 
            if "chr" in line[0]: line[0]=line[0].replace('chr','')
            key = ":".join([line[i] for i in [0,1,3,4]])
            if key not in snplist: continue
            total+=1
            if 'Aloft_Recessive' in snplist[key] or 'Aloft_Dominant' in snplist[key]: lof+=1
            indexGT = line[8].split(':').index("GT")
            for group in snplist[key]:
                for sample in samples:
                    header_idx = header.index(sample)
                    GT = sum([int(i) for i in re.split('[/|]',line[header_idx].split(':')[indexGT]) if i != "."])
                    sample_idx = samples.index(sample)
                    if GT == 1:
                        groups[group][2*sample_idx]+=1
                        CADD_weighted_groups[group][2*sample_idx]+=score[key]
                    if GT == 2:
                        groups[group][2*sample_idx+1]+=1
                        CADD_weighted_groups[group][2*sample_idx+1]+=score[key]
            for gene in genes:
                if line[0] != genes[gene][0]:continue
                if int(line[1]) < genes[gene][1] or int(line[1]) > genes[gene][2]:continue
                total_gene[gene]+=1
                for sample in samples:
                    header_idx = header.index(sample)
                    GT = sum([int(i) for i in re.split('[/|]',line[header_idx].split(':')[indexGT]) if i != "."])
                    sample_idx = samples.index(sample)
                    if GT == 1:
                        kegg_gene[gene][2*sample_idx]+=1
                    if GT == 2:
                        kegg_gene[gene][2*sample_idx+1]+=1
            for pathway in pathways:
                for gene in pathways[pathway]:
                    if line[0] != genes[gene][0]:continue
                    if int(line[1]) < genes[gene][1] or int(line[1]) > genes[gene][2]:continue
                    total_pathway[pathway]+=1
                    for sample in samples:
                        header_idx = header.index(sample)
                        GT = sum([int(i) for i in re.split('[/|]',line[header_idx].split(':')[indexGT]) if i != "."])
                        sample_idx = samples.index(sample)
                        if GT == 1:
                            kegg_pathway[pathway][2*sample_idx]+=1
                        if GT == 2:
                            kegg_pathway[pathway][2*sample_idx+1]+=1
                    break
    

    with open(prefix+'.whole.burden.txt', 'w') as fo:
        fo.write('##Total number of deleterious SNVs:'+'\t'+str(total)+'\n')
        fo.write('##Total number of loss of function SNVs:'+'\t'+str(lof)+'\n')
        fo.write('#sample'+'\t'+'geo'+'\t'+'\t'.join([f"{group}_{stat}" for group in groups for stat in ['Het','Hom']])+'\t'+'\t'.join([f"CADD_weigted_{group}_{stat}" for group in CADD_weighted_groups for stat in ['Het','Hom']])+'\n')
        for sample in samples:
            fo.write(sample+'\t'+prefix+'\t')
            sample_idx = samples.index(sample)
            line=[]
            for group in groups:
                line.append(str(groups[group][2*sample_idx]))
                line.append(str(groups[group][2*sample_idx+1]))
            for group in CADD_weighted_groups:
                line.append(str(CADD_weighted_groups[group][2*sample_idx]))
                line.append(str(CADD_weighted_groups[group][2*sample_idx+1]))
            fo.write('\t'.join(line))
            fo.write('\n')
    
    with open(prefix+'.gene.burden.txt', 'w') as fo:
        for gene in genes:
            fo.write('##Total number of deleterious SNVs of '+gene+':'+'\t'+str(total_gene[gene])+'\n')
        fo.write('#sample'+'\t'+'geo'+'\t'+'\t'.join([f"{gene}_{stat}" for gene in kegg_gene for stat in ['Het','Hom']])+'\n')
        for sample in samples:
            fo.write(sample+'\t'+prefix+'\t')
            sample_idx = samples.index(sample)
            line=[]
            for gene in kegg_gene:
                line.append(str(kegg_gene[gene][2*sample_idx]))
                line.append(str(kegg_gene[gene][2*sample_idx+1]))
            fo.write('\t'.join(line))
            fo.write('\n')

    with open(prefix+'.pathway.burden.txt', 'w') as fo:
        for pathway in pathways:
            fo.write('##Total number of deleterious SNVs of '+pathway+':'+'\t'+str(total_pathway[pathway])+'\n')
        fo.write('#sample'+'\t'+'geo'+'\t'+'\t'.join([f"{pathway}_{stat}" for pathway in kegg_pathway for stat in ['Het','Hom']])+'\n')
        for sample in samples:
            fo.write(sample+'\t'+prefix+'\t')
            sample_idx = samples.index(sample)
            line=[]
            for pathway in kegg_pathway:
                line.append(str(kegg_pathway[pathway][2*sample_idx]))
                line.append(str(kegg_pathway[pathway][2*sample_idx+1]))
            fo.write('\t'.join(line))
            fo.write('\n')

if __name__ == '__main__':
    vcf = sys.argv[1]
    damage = sys.argv[2]
    kegg = sys.argv[3]
    prefix = sys.argv[4]
    BurdenTest(vcf,damage,kegg,prefix) 