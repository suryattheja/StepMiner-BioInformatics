
# coding: utf-8

# In[ ]:


# In[77]:


def stepminer(prefix):
    try:
        expr_df = pd.read_table(prefix+"-expr.txt")
        #print(expr_df)
    except:
        print('-expr file not found')
        exit()
    try:
        idx_df = pd.read_table(prefix+"-idx.txt")
    except:
        print('-idx file not found')
    try:
        ih_df = pd.read_table(prefix+"-ih.txt")
    except:
        print('-ih file not found')
    try:
        surv_df = pd.read_table(prefix+"-survival.txt")
    except:
        print('-survival file not found')
    thresholds=findThresholds(expr_df)
    writeThrTxt(thresholds)
    writeInfoTxt(thresholds)
    writeBVTxt(thresholds)
    
    
    return
        


# In[71]:


def findThresholds(expr_df):
    sample_columns = list(expr_df.columns)[2:]
    thresholds  = []
    for index, row in expr_df.iterrows():
        probe_id = row["ProbeID"]
        name = row["Name"]
        values = row[sample_columns].tolist()
        #print(values)
        thr_vals = fit_step_miner(values)
        thr_vals["name"]=name
        thr_vals["probe_id"]=probe_id
        stat = f_statistic(thr_vals["sse"],thr_vals["sstot"],thr_vals["n"])
        thr_vals["stat"]=round(stat,3)
        thr_vals["thr-0.5"]=round(thr_vals["thr"]-0.5,3)
        thr_vals["thr+0.5"]=round(thr_vals["thr"]+0.5,3)
        thresholds.append(thr_vals)#thr_vals is a dict
        
    return thresholds



# In[72]:


def fit_step_miner(values):  # l is a list
    l=sorted(values)
    # Sort the list if not sorted
    n = len(l)
    n1 = 0
    n2 = n
    min_m1 = m1 = 0
    m2 = sum(l)/n
    min_m2 = m = m2
    min_i = -1
    sse = sum([(x - m2)**2 for x in l])
    sstot = sse
    min_sse = sse
    ssr = 0
    thr = - float('inf')
    for i,x in enumerate(l):
        m1 = ((m1* n1) + x)/(n1+1)
        n1+=1
        if n1 == n:
            m2 = 0
            n2 = 0
        else:
            m2 = ((m2 * n2) - x)/(n2 - 1)
            n2 = n2 - 1
        ssr = (n1 * ((m1 - m)**2)) + (n2 * ((m2 - m)**2))
        sse = (sstot - ssr)
        if min_sse > sse:
            #print(i)
            min_sse = sse
            min_m1 = m1
            min_m2 = m2
            thr = x
            min_i = i
    ret_vals = {}
    #l_mean=np.mean(l)
    l_mean = m
    l_thr=(min_m1+min_m2)/2
    ret_vals["mean"]=round(l_mean,3)
    ret_vals["min"]=round(min(l),3)
    ret_vals["max"]=round(max(l),3)
    ret_vals["mean-thr"]=round(l_mean-l_thr,3)
    ret_vals["sd"]=round(np.std(l),3)
    ret_vals["thr"] =round(l_thr,3)
    ret_vals["sse"] = round(min_sse,3)
    ret_vals["sstot"] = round(sstot,3)
    ret_vals["n"] = n
    
    upperlim=l_thr+0.5
    lowerlim=l_thr-0.5
    BVlist=[]
    for val in values:
        if(val>upperlim):
            BVlist.append('2')
        else: 
            if(val<lowerlim):
                BVlist.append('0')
            else:
                BVlist.append('1')
    BVstring=''.join(BVlist)
    ret_vals["BV"]=BVstring
    
    hi=sum([1 for i in BVlist if i=='2'])
    low=sum([1 for i in BVlist if i=='0'])
    inter=len(BVlist)-hi-low
    #print(values)
    #print(l_thr)
    #print(BVstring)
    ret_vals["hi"]=hi
    ret_vals["lo"]=low
    ret_vals["int"]=inter
    
    #perc
    if(l_thr>l_mean):
        numb_betw=sum([1 for i in l if i>=l_mean and i<=l_thr])
        numb_below=sum([1 for i in l if i<=l_thr])
        perc=-1*numb_betw/numb_below
    else:
        numb_betw=sum([1 for i in l if i<=l_mean and i>=l_thr])
        numb_above=sum([1 for i in l if i>=l_thr])
        perc=numb_betw/numb_above
    ret_vals["perc"]=round(perc,3)
    
    #thrnum
    thrnum=sum([1 for i in l if i<=l_thr])
    ret_vals["thrNum"]=thrnum

    return ret_vals


# In[73]:


def f_statistic(sse,sstot,n):
    ssr = sstot - sse
    dof_ssr = 3 if n > 4 else 2
    dof_sse = n - 4 if n > 4 else 1
    msr = ssr/dof_ssr
    mse = sse/dof_sse
    return(msr/mse)


# In[74]:


def writeThrTxt(thresholds):
    with open(prefix+'-thr.txt','w') as f:
        for dic in thresholds:
                f.write(dic["probe_id"]+"\t"+str(dic["thr"])+"\t"+str(dic["stat"])+"\t"+str(dic["thr-0.5"])+"\t"
                        +str(dic["thr+0.5"])+"\t")
                f.write('\n')
    return


# In[75]:


def writeInfoTxt(thresholds):
    with open(prefix+'-info.txt','w') as f:
        f.write("AffyID"+"\t"+"name"+"\t"+"thr"+"\t"+"mean"+"\t"+"mean-thr"+
                "\t"+"perc"+"\t"+"min"+"\t"+"max"+"\t"+"sd"+"\t"+"thrNum"+"\t"+"hi"+"\t"+"int"+"\t"+"lo"+"\n")
        for t in thresholds:
            f.write(t["probe_id"]+"\t"+t["name"].split(":")[0]+"\t"+str(t["thr"])+"\t"+str(t["mean"])+"\t"+str(t["mean-thr"])+
                "\t"+str(t["perc"])+"\t"+str(t["min"])+"\t"+str(t["max"])+"\t"+str(t["sd"])+"\t"+str(t["thrNum"])+"\t"+str(t["hi"])+"\t"+
                    str(t["int"])+"\t"+str(t["lo"]))
            f.write('\n')
    return


# In[76]:


def writeBVTxt(thresholds):
    with open(prefix+'-bv.txt','w') as f:
        f.write("ArrayID"+"\t"+"name"+"\t"+"BitVector"+"\n")
        for t in thresholds:
            f.write(t["probe_id"]+"\t"+t["name"]+"\t"+t["BV"])
            f.write('\n')
    return
    

