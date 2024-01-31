import pandas as pd  
import numpy as np  
  
# Set constants  
ncutoff = 20  #n carriers you want per window
genename = 'exampleGene'  
  
# Load data  
dfall = pd.read_csv('example_gene_counts', delim_whitespace=True)   #processed from plink hardy command
dfall[['chr', 'pos', 'ref', 'alt']] = dfall['ID'].str.split(':', expand=True)  
  
# Calculate carriers
dfall['carriers'] = np.where(  
    dfall.HOM_A1_CT.astype(int) > dfall.TWO_AX_CT.astype(int),  
    dfall.TWO_AX_CT.astype(int) * 2 + dfall.HET_A1_CT.astype(int),  
    dfall.HET_A1_CT.astype(int) + dfall.HOM_A1_CT.astype(int) * 2  
)  
  
# Separate high frequency variants 
# dfsall has single variant windows
dfsall = dfall[dfall.carriers >= ncutoff]  
dfsall = dfsall.sort_values(by=['chr', 'pos']).reset_index(drop=True)  
  
# Separate variants with multiple carriers at the same site  
# dfmall has windows with only 1 position
dfall = dfall[dfall.carriers < ncutoff]  
dfall['totalcaratsite'] = dfall.groupby(['chr', 'pos']).carriers.transform(np.sum)  
dfmall = dfall[dfall.totalcaratsite >= ncutoff]  
dfmall = dfmall.sort_values(by=['chr', 'pos']).reset_index(drop=True)  
  
# Keep remaining variants in dfall  
dfall = dfall[dfall.totalcaratsite < ncutoff]  
dfall = dfall.sort_values(by=['chr', 'pos']).reset_index(drop=True)  
  
# Filter for gene of interest  
df = dfall[dfall.MANECT_gene == genename].reset_index(drop=True)  
dfs = dfsall[dfsall.MANECT_gene == genename].reset_index(drop=True)  
dfm = dfmall[dfmall.MANECT_gene == genename].reset_index(drop=True)  
  
# Initialize window dataframe  
dfwindows = pd.DataFrame(columns = ['SNP', 'Window', 'carriers', 'solo_variant_window', 'pos', 'solo_site_window'])  
  
# Initialize window counter  
window = 0  
  
# Process windows with only one variant  
for i in range(len(dfs)):  
    dfwindows = dfwindows.append({  
        'SNP': dfs.ID[i],   
        'Window': f'{genename}_{ncutoff}_Window__{window}',  
        'carriers': dfs.carriers[i],  
        'solo_variant_window': 1,  
        'pos': dfs.pos[i],  
        'solo_site_window': 0  
    }, ignore_index=True)  
    window += 1  
  
# Process windows with only one position  
for i in range(len(dfm)):  
    dfwindows = dfwindows.append({  
        'SNP': dfm.ID[i],   
        'Window': f'{genename}_{ncutoff}_Window__{window}',  
        'carriers': dfm.carriers[i],  
        'solo_variant_window': 0,  
        'pos': dfm.pos[i],  
        'solo_site_window': 1  
    }, ignore_index=True)  
    window += 1 if i == len(dfm) - 1 or dfm.pos[i+1] != dfm.pos[i] else 0  
  
# Process windows with multiple sites  
carriers = 0  
# build 1 window per variant, starting at that variant and building until you hit ncutoff carriers
for i in range(len(df)):  
    if carriers != 0:  
        break  
      
    # Do not create a new window if this var is at the same site as the last one
    # (variants at the same site will always be in the same window)
    if i != 0 and df.pos[i] <= df.pos[i - 1]:  
        continue  
  
    # Append current variant to windows dataframe  
    dfwindows = dfwindows.append({  
        'SNP': df.ID[i],  
        'Window': f'{genename}_{ncutoff}_Window__{window}',  
        'carriers': df.carriers[i],  
        'solo_variant_window': 0,  
        'pos': df.pos[i],  
        'solo_site_window': 0  
    }, ignore_index=True)  
  
    # Update total carriers for current window (n carriers in starting variant) 
    carriers += df.carriers[i]  
  
    # Walk through remaining variants in gene until you hit ncutoff carriers
    for x in range(i + 1, len(df)):  
        # Append current variant to windows dataframe  
        dfwindows = dfwindows.append({  
            'SNP': df.ID[x],  
            'Window': f'{genename}_{ncutoff}_Window__{window}',  
            'carriers': df.carriers[x],  
            'solo_variant_window': 0,  
            'pos': df.pos[x],  
            'solo_site_window': 0  
        }, ignore_index=True)  
  
        # Update total carriers for current window with new carriers
        carriers += df.carriers[x]  
  
        # If total carriers for current window reached cutoff, reset carriers and increment window  
        if carriers >= ncutoff:  
            # If the position is the same at the next site, keep adding to this window until
            # there's a new pos
            if x != len(df) - 1 and df.pos[x] == df.pos[x + 1]:  
                continue  
            # If the remaining variants in the rest of the gene have zero carriers,
            # keep adding to this window
            elif np.sum(df[df.index > x].carriers) == 0:  
                continue  
            # If the remaining carriers in the rest of the gene are less than ncutoff,
            # keep adding to this window  
            elif np.sum(df[df.index > x].carriers) + carriers < ncutoff:  
                carriers = 0  
            else:  
                carriers = 0  
                window += 1  
            break
  
# Calculate total carriers per window and sort windows dataframe by position  
dfwindows['carriers_in_window'] = dfwindows.groupby('Window').carriers.transform(np.sum)  
dfwindows = dfwindows.sort_values(by='pos').reset_index(drop=True)  
  
# Write windows dataframe to file  
dfwindows.to_csv('exampleGene_PW.txt', sep="\t", index=False)  
