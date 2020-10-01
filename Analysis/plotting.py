import numpy as np
import math
import freud
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas
import os.path as path
import scipy
from collections import Counter
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl

def create_bins(lower_bound, width, quantity):
    """ create_bins returns an equal-width (distance) partitioning. 
        It returns an ascending list of tuples, representing the intervals.
        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
        and i < quantity, satisfies the following conditions:
            (1) bins[i][0] + width == bins[i][1]
            (2) bins[i-1][0] + width == bins[i][0] and
                bins[i-1][1] + width == bins[i][1]
    """
    

    bins = []
    #edited out a +1 in the second input of the np.arange
    for low in np.arange(lower_bound, 
                     lower_bound + quantity*width, width):
        bins.append((low, low+width))
    return bins

def find_bin(value, bins):
    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        binning returns the smallest index i of bins so that
        bin[i][0] <= value < bin[i][1]
    """
    
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1

def HeatMap(arr1,name1,arr2,name2,constant,namec,df_single_np):
    X, Y = np.meshgrid(arr1, arr2)

    Z = np.zeros_like(Y, dtype=np.float64)


    for i in range(Y.shape[0]):
        for j in range(Y.shape[1]):
            val = df_single_np[(df_single_np[name1]==X[i,j]) & (df_single_np[name2]==Y[i,j])
                        & (np.round(df_single_np[namec],2)==constant)]['fSASA']
            if len(val):
                Z[i,j] = val.values[0]
    return Z

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def Plot_Heatmaps(chainlength,chain_density,fsa,df_single_np,df,cmap,inter):
    markersize=100
    figs, axes = plt.subplots(1,4,figsize=(13,4),gridspec_kw={"width_ratios":[1, 1,1,0.05]})
    #figs.set_clip_on(False)
    #axes[0].set_clip_on(False)
    axes[0].set_autoscale_on(True)
    axes[1].set_autoscale_on(True)
    axes[2].set_autoscale_on(True)

    df_range =  df.loc[df['fSASA']>0]
    norm = mpl.colors.Normalize(vmin=np.min(df_range['fSASA'].values), 
                                vmax=np.max(df_range['fSASA'].values), clip=True)
    
    index=0
    axes[index].set_xlabel('Fractional Surface Area')
    axes[index].set_ylabel('Chain Density, chains / nm^2')
    Z = HeatMap(fsa,'fsa',chain_density,'chain_density',6,'chainlength',df_single_np)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                                 aspect="auto",alpha = 1,
                                 extent=[min(fsa), max(fsa)+0.03, min(chain_density), max(chain_density)], 
                                 norm=norm)
    #,vmin=np.min(sasa),vmax=np.max(sasa))

    booleans_dis=[]
    booleans_str=[]
    booleans_agg=[]
    booleans=[]
    fixed_cl=6
    for cl in df['Chain Length'].values:
        if cl==fixed_cl:
            booleans.append(True)
        else:
            booleans.append(False)
    for phase in df.Phase:
        if phase==0:
            booleans_dis.append(True)
            booleans_str.append(False)
            booleans_agg.append(False)
        if phase==1:
            booleans_dis.append(False)
            booleans_str.append(True)
            booleans_agg.append(False)
        if phase==2:
            booleans_dis.append(False)
            booleans_str.append(False)
            booleans_agg.append(True)
    cl18 = pandas.Series(booleans)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['FSA'],df[stringy]['Chain Density'],
                        marker='o',color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['FSA'],df[agg]['Chain Density'],marker='v',color='red',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['FSA'],df[dis]['Chain Density'],marker='s',color='k',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('Chain Length='+ str(fixed_cl) + ' beads')
    axes[index].set_xlim(0.19,0.66)
    axes[index].set_ylim(2.45,4.55)

    index=1
    axes[index].set_xlabel('Fractional surface area')
    axes[index].set_ylabel('Chain length')
    Z = HeatMap(fsa,'fsa',chainlength,'chainlength',3.5,'chain_density',df_single_np)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                        extent=[min(fsa), max(fsa)+0.03, min(chainlength)-1,
                        max(chainlength)+1], aspect="auto",alpha = 1,norm=norm)
    booleans=[]
    fixed_cd=3.5
    for cd in df['Chain Density'].values:
        if cd==fixed_cd:
            booleans.append(True)
        else:
            booleans.append(False)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['FSA'],df[stringy]['Chain Length'],marker='o',
                        color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['FSA'],df[agg]['Chain Length'],marker='v',color='red',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['FSA'],df[dis]['Chain Length'],marker='s',color='k',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('Chain Density='+ str(fixed_cd) + ' Chains / nm^2')
    axes[index].set_xlim(0.19,0.66)
    axes[index].set_ylim(4.67,11.33)


    index=2
    axes[index].set_xlabel('Chain density, chains / nm^2')
    axes[index].set_ylabel('Chain length')
    Z = HeatMap(chain_density,'chain_density',chainlength,'chainlength',0.55,'fsa',df_single_np)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                        extent=[min(chain_density), max(chain_density), min(chainlength)-1,
                        max(chainlength)+1], aspect="auto",alpha = 1,norm=norm)
    booleans=[]
    fixed_fsa=0.55
    for frac in df['FSA'].values:
        if frac==fixed_fsa:
            booleans.append(True)
        else:
            booleans.append(False)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['Chain Density'],df[stringy]['Chain Length'],marker='o',
                        color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['Chain Density'],df[agg]['Chain Length'],marker='v',color='red',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['Chain Density'],df[dis]['Chain Length'],marker='s',color='k',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('FSA='+ str(fixed_fsa))
    axes[index].set_xlim(2.45,4.55)
    axes[index].set_ylim(4.67,11.33)

    cb = figs.colorbar(heatmap,cax=axes[3])
    cb.set_label('Fractional Solvent-accessible surface area')
    plt.tight_layout()
    plt.savefig('../Figures/heatmap-inter-SASA.pdf',dpi=600)
    
def Plot_Fixed_Fsa(df,fixed_fsa):
    markersize=100
    booleans=[]
    for fsa in df.FSA:
        if fsa==fixed_fsa:
            booleans.append(True)
        else:
            booleans.append(False)
    dis,stringy,agg = Determine_Phase(df,booleans)
    fig,axes=plt.subplots(figsize=(4,4))
    axes.set_xbound(2,5)
    axes.set_ybound(2,5)
    axes.set_ylim(4.67,11.33)
    axes.set_xlim(2.45,4.55)
    axes.scatter(df[stringy]['Chain Density'],df[stringy]['Chain Length'],marker='o',
                 s=markersize,edgecolors = 'k',color='purple',linewidth=1)
    axes.scatter(df[agg]['Chain Density'],df[agg]['Chain Length'],marker='v',s=markersize,edgecolors = 'k',color='red',linewidth=1)
    axes.scatter(df[dis]['Chain Density'],df[dis]['Chain Length'],marker='s',s=markersize,edgecolors = 'k',color='k',linewidth=1)
    axes.set_title('FSA='+ str(fixed_fsa))
    axes.set_xlabel('Chain Density, chains/nm')
    axes.set_ylabel('Chain Length, # beads')
    axes.legend()
    plt.savefig('../Figures/fix_fsa.pdf')
    
def Plot_Fixed_Cd(df,fixed_cd):
    markersize=100
    booleans=[]
    for chain_density in df['Chain Density'].values:
        if chain_density==fixed_cd:
            booleans.append(True)
        else:
            booleans.append(False)
    dis,stringy,agg = Determine_Phase(df,booleans)
    fig,axes=plt.subplots(figsize=(4,4))
    axes.set_xbound(2,5)
    axes.set_ybound(2,5)
    axes.set_xlim(0.19,0.66)
    axes.set_ylim(4.67,11.33)
    axes.scatter(df[stringy]['FSA'],df[stringy]['Chain Length'],marker='o',
                 s=markersize,edgecolors = 'k',color='purple',linewidth=1)
    axes.scatter(df[agg]['FSA'],df[agg]['Chain Length'],marker='v',s=markersize,edgecolors = 'k',color='red',linewidth=1)
    axes.scatter(df[dis]['FSA'],df[dis]['Chain Length'],marker='s',s=markersize,edgecolors = 'k',color='k',linewidth=1)
    axes.set_title('Chain Density='+ str(fixed_cd))
    axes.set_xlabel('FSA')
    axes.set_ylabel('Chain Length, # beads')
    plt.savefig('../Figures/fix_cd.pdf')
    
def Plot_Fixed_Cl(df,fixed_cl):
    markersize=100
    booleans=[]
    for chain_length in df['Chain Length'].values:
        if chain_length==fixed_cl:
            booleans.append(True)
        else:
            booleans.append(False)
    dis,stringy,agg = Determine_Phase(df,booleans)     
    cl18 = pandas.Series(booleans)

    fig,axes=plt.subplots(figsize=(4,4))
    axes.set_xbound(2,5)
    axes.set_ybound(2,5)
    axes.set_xlim(0.2,0.66)
    axes.set_ylim(2.45,4.55)
    axes.scatter(df[stringy]['FSA'],df[stringy]['Chain Density'],marker='o',
                 s=markersize,edgecolors = 'k',color='purple',linewidth=1)
    axes.scatter(df[agg]['FSA'],df[agg]['Chain Density'],marker='v',s=markersize,edgecolors = 'k',color='red',linewidth=1)
    axes.scatter(df[dis]['FSA'],df[dis]['Chain Density'],marker='s',s=markersize,edgecolors = 'k',color='k',linewidth=1)
    axes.set_title('Chain Length='+ str(fixed_cl))
    axes.set_xlabel('FSA')
    axes.set_ylabel('Chain Density, chains/nm')
    plt.savefig('../Figures/fix_cl.pdf')
    
def Determine_Phase(df,booleans):
    booleans_dis=[]
    booleans_str=[]
    booleans_agg=[]
    for phase in df.Phase:
        if phase==0:
            booleans_dis.append(True)
            booleans_str.append(False)
            booleans_agg.append(False)
        if phase==1:
            booleans_dis.append(False)
            booleans_str.append(True)
            booleans_agg.append(False)
        if phase==2:
            booleans_dis.append(False)
            booleans_str.append(False)
            booleans_agg.append(True)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    return(dis,stringy,agg)

def Total_Counts(df,bins,attribute):
    bin_center = []
    for bin in bins:
        center = (bin[0]+bin[1])/2
        bin_center.append(center)

    binned_sasa = []
    for value in df[(df[attribute]>0)][attribute].values:
        bin_index = find_bin(value, bins)
        binned_sasa.append(bin_index)

    frequencies = Counter(binned_sasa)

    counts = []
    for index, center in enumerate(bin_center):
        counts.append(frequencies[index])
    return(counts)

def Plot_Phase(df,phase,total_counts,bins,attribute,in_place=True):
    bin_center = []
    for bin in bins:
        center = (bin[0]+bin[1])/2
        bin_center.append(center)
    binned_sasa = []
    for index,value in enumerate(df[(df[attribute]>0)][attribute].values):
        if df[(df[attribute]>0)]['Phase'].values[index] == phase:
            bin_index = find_bin(value, bins)
            binned_sasa.append(bin_index)

    frequencies = Counter(binned_sasa)
    counts = []
    for index, center in enumerate(bin_center):
        counts.append(frequencies[index])
    for index,value in enumerate(counts):
        if value:
            counts[index] = (value/total_counts[index])
    if in_place:
        plt.plot(bin_center,counts,'b')
    return(counts)

def Calc_Coordnum(data,r_cut):
    neighbor_arr=np.zeros((data.shape[0],25))
    frame_num=0
    for frame in (data):
        for i in np.arange(0,25,1):
            neighbor=0
            for j in np.arange(0,25,1):
                if i!=j:
                    array=frame[i,:]-frame[j,:]
                    #print('Array: ',array)
                    distance=np.sqrt(np.dot(array,array))
                    if distance<=r_cut:
                        neighbor+=1 
            neighbor_arr[frame_num,i]+=neighbor
        frame_num+=1 
    return(np.average(neighbor_arr,axis=0))

def plot_rdf(points_arr, prop, box, r_max=12, bins=200, label=None, ax=None):
    """Helper function for plotting RDFs."""
    box=freud.box.Box.from_box([101, 101, 101])
    #if ax is None:
        #fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        #ax.set_title(prop, fontsize=16)
    rdf = freud.density.RDF(bins,r_max,r_min=4)
    for points in points_arr:
        wrapped_points=box.wrap(points)
        rdf.compute(system=(box, wrapped_points), reset=False)
    if label is not None:
        #ax.plot(rdf.bin_centers, getattr(rdf, prop), label=label)
        ax.legend()
    else:
        counts=getattr(rdf, prop)
        #ax.plot(rdf.bin_centers, getattr(rdf, prop)/counts.max())
        counts=getattr(rdf, prop)
        print(counts.max())
    return (ax,rdf,counts,prop)

def Wrap_data(data,box_length):
    wrapped_data=np.zeros((data.shape[0],25,3))
    wrapped_data[:,:,0]=data[:,:,0]+box_length/2
    wrapped_data[:,:,1]=data[:,:,1]+box_length/2
    wrapped_data[:,:,2]=data[:,:,2]+box_length/2
    for frame in np.arange(0,data.shape[0]-1,1):
        for i in np.arange(0,25,1):
            for j in np.arange(0,25,1):
                if i!=j:
                    array=wrapped_data[frame,i,:]-wrapped_data[frame,j,:]
                    distance=np.sqrt(np.dot(array,
                                            array))
    return(wrapped_data)

def Plot_All_Phases(attribute,df,bins): 
    total_counts=Total_Counts(df,bins,attribute)
    bin_center = []
    for bin in bins:
        center = (bin[0]+bin[1])/2
        bin_center.append(center)
    counts0 = Plot_Phase(df,0,total_counts,bins,attribute,in_place=False)
    counts1 = Plot_Phase(df,1,total_counts,bins,attribute,in_place=False)
    counts2 = Plot_Phase(df,2,total_counts,bins,attribute,in_place=False)
    fig,ax=plt.subplots(1,1,figsize=(6,4))
    ax.plot(bin_center,counts0,ls='-',c='k',label='Dispersed')
    ax.plot(bin_center,counts1,ls='--',c='indigo',label='Stringy')
    ax.plot(bin_center,counts2,ls=':',c='red',label='Aggregated')
    ax.fill_between(bin_center, 0, counts0,facecolor='k', interpolate=True,alpha=0.25)
    ax.fill_between(bin_center, 0, counts1,facecolor='indigo', interpolate=True,alpha=0.25)
    ax.fill_between(bin_center, 0, counts2,facecolor='red', interpolate=True,alpha=0.25)
    ax.set_ylim(0,1.1)
    #plt.legend(title='Phase',fontsize=8)
    plt.xlabel(attribute)
    plt.ylabel('Normalized Phase')
    plt.savefig('../Figures/histogram-' + attribute +'.pdf')
    plt.show()
    
def Plot_Heatmaps_all(chainlength,chain_density,fsa,df_single_np,df,name_overlay,cmap,inter):
    markersize=100
    figs, axes = plt.subplots(1,4,figsize=(13,4),gridspec_kw={"width_ratios":[1, 1,1,0.05]})
    #figs.set_clip_on(False)
    #axes[0].set_clip_on(False)
    axes[0].set_autoscale_on(True)
    axes[1].set_autoscale_on(True)
    axes[2].set_autoscale_on(True)
    
    df_range =  df.loc[df['fSASA']>0]
    norm = mpl.colors.Normalize(vmin=np.min(df_range[name_overlay].values), 
                                vmax=np.max(df_range[name_overlay].values), clip=False)
    
    index=0
    axes[index].set_xlabel('Fractional Surface Area')
    axes[index].set_ylabel('Chain Density, chains / nm^2')
    Z = HeatMap_all(fsa,'fsa',chain_density,'chain_density',6,'chainlength',df_single_np,name_overlay)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                                 aspect="auto",alpha = 1,
                                 extent=[min(fsa), max(fsa)+0.03, min(chain_density), max(chain_density)], 
                                 norm=norm)
    #,vmin=np.min(sasa),vmax=np.max(sasa))

    booleans_dis=[]
    booleans_str=[]
    booleans_agg=[]
    booleans=[]
    fixed_cl=6
    for cl in df['Chain Length'].values:
        if cl==fixed_cl:
            booleans.append(True)
        else:
            booleans.append(False)
    for phase in df.Phase:
        if phase==0:
            booleans_dis.append(True)
            booleans_str.append(False)
            booleans_agg.append(False)
        if phase==1:
            booleans_dis.append(False)
            booleans_str.append(True)
            booleans_agg.append(False)
        if phase==2:
            booleans_dis.append(False)
            booleans_str.append(False)
            booleans_agg.append(True)
    cl18 = pandas.Series(booleans)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['FSA'],df[stringy]['Chain Density'],
                        marker='o',color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['FSA'],df[agg]['Chain Density'],marker='v',color='red',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['FSA'],df[dis]['Chain Density'],marker='s',color='k',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('Chain Length='+ str(fixed_cl) + ' beads')
    axes[index].set_xlim(0.19,0.66)
    axes[index].set_ylim(2.45,4.55)

    index=1
    axes[index].set_xlabel('Fractional surface area')
    axes[index].set_ylabel('Chain length')
    Z = HeatMap_all(fsa,'fsa',chainlength,'chainlength',3.5,'chain_density',df_single_np,name_overlay)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                        extent=[min(fsa), max(fsa)+0.03, min(chainlength)-1,
                        max(chainlength)+1], aspect="auto",alpha = 1,norm=norm)
    booleans=[]
    fixed_cd=3.5
    for cd in df['Chain Density'].values:
        if cd==fixed_cd:
            booleans.append(True)
        else:
            booleans.append(False)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['FSA'],df[stringy]['Chain Length'],marker='o',
                        color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['FSA'],df[agg]['Chain Length'],marker='v',color='red',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['FSA'],df[dis]['Chain Length'],marker='s',color='k',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('Chain Density='+ str(fixed_cd) + ' Chains / nm^2')
    axes[index].set_xlim(0.19,0.66)
    axes[index].set_ylim(4.67,11.33)


    index=2
    axes[index].set_xlabel('Chain density, chains / nm^2')
    axes[index].set_ylabel('Chain length')
    Z = HeatMap_all(chain_density,'chain_density',chainlength,'chainlength',0.55,'fsa',df_single_np,name_overlay)
    heatmap = axes[index].imshow(Z, cmap=cmap, interpolation=inter, origin='lower',
                        extent=[min(chain_density), max(chain_density), min(chainlength)-1,
                        max(chainlength)+1], aspect="auto",alpha = 1,norm=norm)
    booleans=[]
    fixed_fsa=0.55
    for frac in df['FSA'].values:
        if frac==fixed_fsa:
            booleans.append(True)
        else:
            booleans.append(False)
    dis = pandas.Series(np.logical_and(booleans,booleans_dis))
    stringy = pandas.Series(np.logical_and(booleans,booleans_str))
    agg = pandas.Series(np.logical_and(booleans,booleans_agg))
    axes[index].scatter(df[stringy]['Chain Density'],df[stringy]['Chain Length'],marker='o',
                        color='purple',edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[agg]['Chain Density'],df[agg]['Chain Length'],marker='v',color='red',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].scatter(df[dis]['Chain Density'],df[dis]['Chain Length'],marker='s',color='k',
                        edgecolors = 'white',s=markersize,linewidth=1.5)
    axes[index].set_title('FSA='+ str(fixed_fsa))
    axes[index].set_xlim(2.45,4.55)
    axes[index].set_ylim(4.67,11.33)


    cb = figs.colorbar(heatmap,cax=axes[3])
    cb.set_label(name_overlay)
    plt.tight_layout()
    plt.savefig('../Figures/heatmap-newmap.pdf',dpi=600)
            
    
def HeatMap_all(arr1,name1,arr2,name2,constant,namec,df_single_np,name_of_interest):
    X, Y = np.meshgrid(arr1, arr2)

    Z = np.zeros_like(Y, dtype=np.float64)


    for i in range(Y.shape[0]):
        for j in range(Y.shape[1]):
            val = df_single_np[(df_single_np[name1]==X[i,j]) & (df_single_np[name2]==Y[i,j])
                        & (np.round(df_single_np[namec],2)==constant)][name_of_interest]
            if len(val):
                Z[i,j] = val.values[0]
    return Z

def Get_Clrmap(c_phases,meshing,bounds):
    regions = len(c_phases)
    colorlist = []
    intermediates = []
    for i,bound in enumerate(bounds[:-1]):
        intermediates.append((bounds[i+1]-bound)/2+bound)
    intermediates = np.array(intermediates)
    scale = 1/(bounds.max()-bounds.min())
    
    #color region 1
    steps = np.floor((bounds[1]-bounds[0])*scale*meshing)
    #inc = (c_phases[1] - c_phases[0])/steps
    for a in np.arange(0,steps):
        c_val = c_phases[0]
        colorlist.append(c_val)
    
    #region intermediate
    steps = np.floor((intermediates[1]-bounds[1])*scale*meshing) 
    inc = (c_phases[1] - c_phases[0])/steps
    for a in np.arange(0,steps):
        c_val = c_phases[0] + a*inc
        colorlist.append(c_val)
        
    steps = np.floor((bounds[2]-intermediates[1])*scale*meshing) 
    inc = (c_phases[2] - c_phases[1])/steps
    for a in np.arange(0,steps):
        c_val = c_phases[1] + a*inc
        colorlist.append(c_val)
        
    #color region 2
    steps = np.floor((bounds[3]-bounds[2])*scale*meshing)
    for a in np.arange(0,steps):
        c_val = c_phases[2]
        colorlist.append(c_val)
        
    #region intermediate
    steps = np.floor((intermediates[3]-bounds[3])*scale*meshing)
    inc = (c_phases[3] - c_phases[2])/steps
    for a in np.arange(0,steps):
        c_val = c_phases[2] + a*inc
        colorlist.append(c_val)
        
    steps = np.floor((bounds[4]-intermediates[3])*scale*meshing) 
    inc = (c_phases[4] - c_phases[3])/steps
    for a in np.arange(0,steps):
        c_val = c_phases[3] + a*inc
        colorlist.append(c_val)
    
    
    #color region 3
    steps = np.floor((bounds[5]-bounds[4])*scale*meshing)
    for a in np.arange(0,steps):
            c_val = c_phases[4]
            colorlist.append(c_val)
    return(colorlist)