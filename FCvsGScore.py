
import numpy as np
import pandas as pan
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_rows(file):
    fin = open(file, "r")
    buf = fin.read()
    fin.close()
    rows = buf.split("\n")
    return rows

def get_data(datafile, xname, yname):
	#loading dataset
	df = pan.read_csv (datafile, sep = ";")
	
	df["region"] = df ["region"].map({"CdsExon" : "red","UtrExon3" : "blue", "UtrExon5" : "green", "Intron" : "yellow", "Downstream":"purple", "Upstream": "cyan"})
	print(df["region"])
	return df

def completeinformation(datafile, infile, outfile):
    
    first = True
    rows = get_rows(datafile)
    
    genes = []
    codes = []
    fcs = []
    pvals = []
    #padjs = []
    
    for i in range(1,len(rows)):
        cols = rows[i].split(";")
        if len(cols) > 1:
            genes.append(cols[0])
            codes.append(cols[0])
            fcs.append(cols[1])
            pvals.append(cols[2])
    #        padjs.append(cols[6])
    
    
    headbuf = ["gene","code","region","#","Position","Length","QGRS","G-Score"]        
    header = ["yvar","pvar"]
    
    regions = ["CdsExon", "Intron", "Upstream", "Downstream", "UtrExon3", "UtrExon5"]
    df = pan.read_csv(infile, sep = ";",index_col=False)
    df = df[["gene","code","region","#","Position","Length","QGRS","G-Score"]]
    df_out = pan.DataFrame()
    
    for i in range(len(codes)):
       
        buf = df[(df.code==codes[i])]   
        count = buf.shape[0]
        listadd = []
        listadd.append(str(fcs[i]))
        listadd.append(str(pvals[i]))
        #listadd.append(str(padjs[i]))
        listindex = []
        for n in range(count):
            listindex.append(n)
        buf = pan.DataFrame(data = buf.values, index = listindex, columns = headbuf)
        df_const = pan.DataFrame(count*[listadd], columns = header)
        print(buf)
        print (df_const)
        df_buf = pan.DataFrame.join(buf, df_const)
        if first == True:
            df_buf.to_csv(outfile, index = False, header = True,sep =  ";")
            first = False
        else:
            df_buf.to_csv(outfile,mode = "a", index = False, header = False,sep =  ";")
    return 
    
    rows = qlibs.get_rows(datafile)
    
    codes = []
    yvars = []
    pvar = []
    
    
    for row in rows:
        cols = row.split(";")
        if len(cols) > 1:
            #genes.append(cols[0])
            codes.append(cols[0])
            yvars.append(cols[1])
            pvars.append(cols[2])
            
    
    text = "gene;code;region;#;Position;Length;QGRS;G-Score;;log2FoldChange;pvalue;padj;\n"
    #text = "gene;code;region;#;"";"start";"width";"strand";"score";"seq";"
    
    rows = qlibs.get_rows(infile)
    for row in rows:
        cols = row.split(";")
        if len(cols) > 1:
            k = qlibs.find(cols[1], codes)
            if k >= 0:
                text = text + row.replace(";;;", ";;") + fcs[k] + ";" + pvals[k] + ";" + padjs[k] + ";\n"
    
    fout = open(outfile,"w")
    fout.write(text)
    fout.close()
    
    return 

def colorgraph(datafile, xname, yname, title, no, show = False):	
    df = get_data(datafile, xname, yname)
    x = df[xname]
    y = df[yname]
    y = y.apply(lambda x: float(x.replace(',','.')))

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set(xlabel=xname, ylabel=yname)
    ax.scatter(x,y,alpha = 0.70, c = df['region'], cmap = cm.brg)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color = "grey", linestyle ="-", linewidth = 0.25, alpha =0.5)
            
    patch1 = mpatches.Patch(color='red', label='CdsExon')	
    patch2 = mpatches.Patch(color='blue', label='UtrExon3')
    patch3 = mpatches.Patch(color='green', label='UtrExon5')
    patch4 = mpatches.Patch(color='yellow', label='Intron')
    patch5 = mpatches.Patch(color='purple', label='Downstream')
    patch6 = mpatches.Patch(color='cyan', label='Upstream')

  

    drawname = datafile.replace(".csv",".png")
    plt.savefig(drawname)
    plt.show()
    
    if show == True: plt.show()
    return    

def pointgraph(datafile, xname, yname, title, no, show = False):	
    df = get_data(datafile, xname, yname)
    x = df[xname]
    y = df[yname]
    y = y.apply(lambda x: float(x.replace(',','.')))

    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set(xlabel=xname, ylabel=yname)
    ax.scatter(x,y,alpha = 0.70)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color = "grey", linestyle ="-", linewidth = 0.25, alpha =0.5)
            

    
    drawname = datafile.replace(".csv",".png")
    plt.savefig(drawname)
    plt.show()
    
    if show == True: plt.show()
    return   


###################
datafile = "C:\\Users\\Utente\\Desktop\\FCvsGScore\\datafile.csv"
infile = "C:\\Users\\Utente\\Desktop\\FCvsGScore\\infile.csv"
outfile = "C:\\Users\\Utente\\Desktop\\FCvsGScore\\outfile.csv"

completeinformation(datafile, infile, outfile)

#############
score = "G-Score"
value= 'yvar' #'log2FoldChange'
title= 'A' 
file = "C:\\Users\\Utente\\Desktop\\FCvsGScore\\outfile.csv"
pointgraph(file, score, value, title, True)
colorgraph(file, score, value, title, True)
############
