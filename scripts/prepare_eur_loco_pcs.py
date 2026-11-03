from pathlib import Path
import os,subprocess,json,hashlib
import numpy as np,pandas as pd
import argparse
parser=argparse.ArgumentParser(description='Compute five genome-wide leave-one-chromosome-out genotype PCs for the expression-matched LDREF cohort.')
parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
root=parser.parse_args().project_root.resolve()
work=root/'results/pca_preparation';work.mkdir(parents=True,exist_ok=True)
out=root/'data/covariates/eur_loco';out.mkdir(parents=True,exist_ok=True)
plink=root/'tools'/('plink2.exe' if os.name == 'nt' else 'plink2')
cols=pd.read_csv(root/'data/expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt',sep='\t',nrows=0).columns
fam=pd.read_csv(root/'data/LDREF/1000G.EUR.1.fam',sep=r'\s+',header=None,dtype=str)
samples=fam.loc[fam[1].isin(cols),[0,1]].copy();samples.columns=['FID','IID'];assert len(samples)==344
keep=work/'keep.txt';samples.to_csv(keep,sep='\t',header=False,index=False)
matrices={};counts={};commands=[];markers=[]
def run(args):
 cmd=[str(plink),*map(str,args)];commands.append(cmd)
 r=subprocess.run(cmd,capture_output=True,text=True)
 if r.returncode:raise RuntimeError(r.stdout+r.stderr)
for c in range(1,23):
 prefix=root/f'data/LDREF/1000G.EUR.{c}';p=work/f'chr{c}'
 run(['--bfile',prefix,'--keep',keep,'--snps-only','just-acgt','--max-alleles','2','--maf','.05','--geno','0','--rm-dup','exclude-all','--indep-pairwise','200','50','.2','--out',p])
 run(['--bfile',prefix,'--keep',keep,'--extract',str(p)+'.prune.in','--export','A','--out',p])
 raw=pd.read_csv(str(p)+'.raw',sep=r'\s+',dtype={'FID':str,'IID':str}).set_index('IID').loc[samples.IID]
 assert raw.FID.tolist()==samples.FID.tolist()
 X=raw.iloc[:,5:].to_numpy(dtype=float);assert np.isfinite(X).all()
 f=X.mean(axis=0)/2;Z=(X-2*f)/np.sqrt(2*f*(1-f))
 matrices[c]=Z@Z.T;counts[c]=Z.shape[1]
 markers.extend({'chrom':c,'snp':s} for s in Path(str(p)+'.prune.in').read_text().splitlines())
 print('PCA markers prepared:',c,Z.shape[1],flush=True)
total=sum(matrices.values());N_MARKERS=sum(counts.values())
for c in range(1,23):
 K=(total-matrices[c])/(N_MARKERS-counts[c]);values,vectors=np.linalg.eigh(K)
 pcs=vectors[:,np.argsort(values)[::-1][:5]]
 df=samples.copy()
 for i in range(5):df[f'PC{i+1}']=pcs[:,i]
 df.to_csv(out/f'eur_loco_chr{c}.tsv',sep='\t',index=False)
pd.DataFrame(markers).to_csv(out/'pca_markers.tsv.gz',sep='\t',index=False,compression='gzip')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
manifest={'method':'genotype PCA on frequency-standardized LD-pruned markers; leave tested chromosome out',
 'samples':344,'PCs':5,'markers_total':N_MARKERS,'markers_by_chromosome':counts,
 'marker_QC':'MAF >= .05; no missing genotypes; biallelic A/C/G/T SNPs; duplicate IDs removed',
 'LD_pruning':'200 variants, 50 variant step, r2 .2','association_variants':'unchanged; marker filters apply only to PCA',
 'commands':commands,'genotype_sha256':{str(p.relative_to(root)):sha(p) for p in sorted((root/'data/LDREF').glob('*')) if p.suffix in ['.bed','.bim','.fam']}}
(out/'sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('All 22 LOCO PC tables saved.',N_MARKERS,flush=True)
