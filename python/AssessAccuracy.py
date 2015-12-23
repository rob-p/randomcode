import matplotlib
matplotlib.use('Agg')
import pandas as pd
import click

def loadRSEM(tf):
    t = pd.read_table(tf, sep='\t', header=True)
    d = {}
    for c in t.columns:
        if c == 'TPM':
            d[c] = 'TPM_true'
        else:
            d[c] = c
    t.rename(columns=d, inplace=True)
    t.set_index("Name", inplace=True)
    return t

def loadFlux(tf):
    t = pd.read_table(tf, sep='\t', header=None, names=None)
    print(t.head())
    t['count'] = t[9] / 2
    #t['TPM_true'] = 10**6 * t[5]
    ts = t.loc[t[3] > 0, :]
    print(len(ts))
    t['TPM_true'] = 10**6 * ((ts['count'] / ts[3]) / (ts['count'] / ts[3]).sum())
    t['Name'] = t[1]
    t.set_index("Name", inplace=True)
    return t

def loadTruth(tf):
    h = open(tf).readline()
    if h.split() == ['transcript_id', 'gene_id', 'length', 'effective_length', 'count', 'TPM', 'FPKM', 'IsoPct']:
        return loadRSEM(tf)
    else:
        return loadFlux(tf)

def loadKal(ef):
    p = pd.read_table(ef, sep='\t')
    p['Name'] = p['target_id']
    p['NumReads'] = p['est_counts']
    p['TPM'] = p['tpm']
    p.set_index("Name", inplace=True)
    return p

def loadSF(ef):
    p = pd.read_table(ef, sep='\t')
    p.set_index('Name', inplace=True)
    return p

def loadEst(ef):
    h = open(ef).readline()
    if h.split() == ['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads']:
        return loadSF(ef)
    else:
        return loadKal(ef)

@click.command()
@click.option('--pred', help='file containing predicted abundances')
@click.option('--true', help='file containing the true abundances.')
def main(pred, true):
    print("pred = {}; true={}".format(pred, true))
    p = loadEst(pred)
    t = loadTruth(true)
    m = t.join(p)

    print("Pearson = {}".format(m['count'].corr(m['NumReads'])))
    print("Spearman = {}".format(m['count'].corr(m['NumReads'], method='spearman')))
    print("Pearson TPM = {}".format(m['TPM_true'].corr(m['TPM'])))
    print("Spearman TPM = {}".format(m['TPM_true'].corr(m['TPM'], method='spearman')))

    from matplotlib import pyplot as plt
    import seaborn as sns
    import numpy as np
    sns.jointplot(np.log(m['TPM_true']+1.0), np.log(m['TPM']+1.0))
    plt.savefig('tpm_corr.pdf')

if __name__ == "__main__":
    main()
