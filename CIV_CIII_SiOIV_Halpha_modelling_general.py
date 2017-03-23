# to run python CIV_CIII_SiOIV_modelling_general.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name # remember to after 'path_to_data' to include the symbol '/'
import matplotlib
matplotlib.use('Agg') 
#matplotlib.use('Agg')
from matplotlib import pylab
import astropy.units as un
import numpy as np
import pyspeckit
import os
import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.interpolate import UnivariateSpline as Interp
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from astropy.cosmology import FlatLambdaCDM
from fitcode.line import *
import json
from matplotlib import rcParams
import time
plot_best_fits=1 # 1 if you want to plot the best fits. to run in cluster must be set to 0
treshold=0.05 #this is the percentage of most negative pixels that are going to be exclude in the fit after the continuum removal
              # this procedure is done iteratively for 3 times.
              #For low-z quasars with little absortion treshold=0.03 is ok (z<3). For  high-z quasars treshold=0.05 would be better
use_exclude_files=0 # 0 if you want the code to fit the lines witouth using manually selected continuum windwos and manually selected 
                    #absorption features. The files should be in the folder ./excludes/ with the following filenames:
                    # exclude_cont_CIV.txt and exclude_CIV.txt. Replace CIV for the name of other lines
                    # you want to fit e.g. CIII, SiOIV, Lya. 



rcParams['text.usetex'] = True #Very important to force python to recognize Latex code


#------Executing the configuration file that includes the requiered parameters---#
execfile("./constraints_single.cfg") 
#--- If you want to change width limits, line center limits, add gaussian components, this is the file to configure---#


if len(sys.argv)<4:
    data_dir='./data/spectra/' #dir related with the spectra that will be fit
    fn='./data/filenames.txt' #file with the list of name of the spectra located in data_dir
else:

    data_dir=sys.argv[2]
    fn=sys.argv[3]


filenames=np.genfromtxt(fn,unpack=1,dtype='str')






if len(sys.argv)<2:
    line_to_fit='CIV'
else:
    line_to_fit=sys.argv[1]
    
    
print line_to_fit
dataout_path=data_dir+'fit/'
plots_path= data_dir+'plots/'



xmin=1100
xmax=9000
n_spec=1



if not os.path.exists(dataout_path):
    os.mkdir(dataout_path)
if not os.path.exists(plots_path):
    os.mkdir(plots_path)
#pylab.ioff()

sample_dictionary={}
try:
    if os.path.isfile( dataout_path + "sample_dictionary.json"):    
        with open(dataout_path + 'sample_dictionary.json', 'rb') as fp:
            sample_dictionary=json.load(fp)
except:
    pass
    




if line_to_fit in ['CIV','CIII','SiOIV','Lya','C']:
    spec_division=['FUV']
    if line_to_fit=='CIV': 
        xminc=1450
        xmaxc=1723        
        
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_CIV.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        exclude_cont=[1000,1445,1465,1695,1725,10000]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]
    if line_to_fit=='CIII': 
        xminc=1678
        xmaxc=2017
        exclude_cont=[1000,1695,1725,1960,2020,10000]
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_CIII.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]
    if line_to_fit=='SiOIV': 
        xminc=1320
        xmaxc=1480
        
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_SiOIV.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        exclude_cont=[1000,1340,1360,1430,1460,10000]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]
    if line_to_fit=='Lya': 
        print 'fitting Ly-alpha'
        xminc=1100
        xmaxc=1480
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_Lya.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        exclude_cont=[1000,1270,1280,1420,1460,10000]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]
    if line_to_fit=='C': 
        xminc=1450
        xmaxc=2017
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_C.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        exclude_cont=[1000,1430,1460,1960,2020,10000]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]


if line_to_fit in ['Halpha','Hbeta']:
    spec_division=['OP']
    if line_to_fit=='Halpha':
        print 'Halpha'
        xminc=6000.0
        xmaxc=7000.0
        if use_exclude_files:
            exclude_file = "./excludes/exclude_cont_Halpha.txt"
            exclude_conti=np.loadtxt(exclude_file,skiprows=2)
        exclude_cont=[1000,6198,6215,6880,6920,10000]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]









try:
    sp_to_run=range(len(filenames))
except:
    filenames=np.array([str(filenames)])
    print 'only one file'
    sp_to_run=[0]


#sp_to_run=[0]
for i in sp_to_run:
    t0=time.time()
    

    fileroot= filenames[i]
    spectrum_file=file=data_dir+fileroot
    print 'fitting ',line_to_fit,' in ', fileroot,'\n\n\n\n'
    
    object_dictionary={}
    object_dictionary1={}
    

        
    plot_objpath=plots_path
    if not os.path.exists(plot_objpath):
        os.mkdir(plot_objpath)

    
    

    

    sp = pyspeckit.Spectrum(spectrum_file)

    if use_exclude_files:            
        exclude_cont=exclude_conti[i,:][:]
        xminc1=exclude_cont[1]
        xmaxc1=exclude_cont[-2]
        print xminc1, xmaxc1, 'cont limits!!!!'
    if line_to_fit=='Lya':
        xminc1=1100
	xmaxc1=1460
    argxmin=np.argmin(np.abs(sp.xarr.value-xminc1))
    argxmax=np.argmin(np.abs(sp.xarr.value-xmaxc1))
    try:
        sp.crop(argxmin,argxmax)
    except:
        continue
    continuous=0.0*sp.data 

    CIV_fit=0.0*sp.data 
    CIII_fit=0.0*sp.data 
    C_fit=0.0*sp.data 
    
    SiOIV_fit=0.0*sp.data 
    
    Lya_fit=0.0*sp.data 

    Halpha_fit=0.0*sp.data 
    
    #------------opening spectrum i---------------------#
    
        
    # -----------set up units properly------------------#        
    
    try:
        mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    except:
        continue
    
    sp.xarr.unit==un.Angstrom
    #sp.xarr.units='angstrom'
    sp.xarr.xtype = 'wavelength'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.data = 10**(-1.0*mag_order)*np.array(sp.data.tolist())    
    sp.error= 10**(-1.0*mag_order)*np.array(sp.error.tolist())    
    #-------------- set up unit properly------------#
    copy=sp.copy()
    if copy.error.mean()==0:
        copy.error=copy.data.mean()*0.1*copy.data
    #-------------continuous fitting-------------------#
    wlmin=sp.xarr[0]
    wlmax=sp.xarr[-1]

    
    wlimin_FUV=xminc
            
            
    #exclude_cont=np.loadtxt(exclude_file,skiprows=2)
    #exclude_cont=exclude_cont[i,:][:]

    backup=sp.copy()
    #backup.plotter(xmin=xminc1,xmax=xmaxc1)
    backup.baseline.powerlaw=False
    backup.baseline(xmin=xminc1, xmax=xmaxc1, exclude=exclude_cont, subtract=False, reset_selection=False, highlight_fitregion=False,powerlaw=False,quiet=True,LoudDebug=False,annotate=False)
            
            
    continuous=backup.baseline.basespec
    continuous_FUV=continuous
    

    
    
    cut=np.percentile(continuous/sp.data,99.5)
    ratio=continuous/sp.data
    
    
    
    

    pylab.rcParams["figure.figsize"]=16,8


    


            
        
        
        
        

            
            

            

    

    
    
    #sp.data[:index]=sp.data[:index] - balmer_template.data[:index] 
    #continuous,wlmin_FUV,L_model=continuous_substraction( i, sp, mag_order,FUV_limits,w)
    
        

    #-------------continuous subtraction-------------------#
    argmax=np.argmax(sp.data)
    posmax=sp.xarr.value[argmax]
    fluxmax=sp.data[argmax]

    dv=10000 #km/s maximum negative velocity to look for BAL features
    dv0=0 #km/s miminum negative velocity to look for BAL features

    dv1=20000 
    dv2=2000

    dv3=1000
    dv4=-1000
    
    dv5=3000
    dv6=-3000

    
    dv7=20000
    dv8=10000

    

    dx=dv*posmax/3e5
    dx0=dv0*posmax/3e5
    dx1=dv1*posmax/3e5
    dx2=dv2*posmax/3e5
    dx3=dv3*posmax/3e5
    dx4=dv4*posmax/3e5
    dx5=dv5*posmax/3e5
    dx6=dv6*posmax/3e5
    dx7=dv7*posmax/3e5
    dx8=dv8*posmax/3e5
    

    poslim=posmax-dx
    poslim0=posmax-dx0
    poslim1=posmax-dx1
    poslim2=posmax-dx2
    poslim3=posmax-dx3
    poslim4=posmax-dx4
    poslim5=posmax-dx5
    poslim6=posmax-dx6
    poslim7=posmax-dx7
    poslim8=posmax-dx8
    

    arglim=np.argmin(np.abs(sp.xarr.value-poslim))
    arglim0=np.argmin(np.abs(sp.xarr.value-poslim0))
    arglim1=np.argmin(np.abs(sp.xarr.value-poslim1))
    arglim2=np.argmin(np.abs(sp.xarr.value-poslim2))
    arglim3=np.argmin(np.abs(sp.xarr.value-poslim3))
    arglim4=np.argmin(np.abs(sp.xarr.value-poslim4))
    arglim5=np.argmin(np.abs(sp.xarr.value-poslim5))
    arglim6=np.argmin(np.abs(sp.xarr.value-poslim6))
    arglim7=np.argmin(np.abs(sp.xarr.value-poslim7))
    arglim8=np.argmin(np.abs(sp.xarr.value-poslim8))
    

    sp.data=sp.data - continuous_FUV
    SN=np.median(copy.data/sp.error)


    print 'signal to noise ratio =',SN
    if line_to_fit=='CIV':
        
        limits_CIV[3]=(0,sp.data.max()) #Defining CIV ampplitud limits
        limits_CIV[0]=(0,sp.data.max()/3.0) #Defining CIV ampplitud limits
        limited_CIV[0]=(True,True) #Defining CIV ampplitud limits
        guesses_CIV[3]=sp.data.max()
        guesses_CIV[6]=sp.data.max()/3.0
        guesses_CIV[9]=sp.data.max()
        guesses_CIV[12]=sp.data.max()/3.0
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
            
        #exclude=[1489.2,1495.6,1497.73,1518.9,1521.0,1532.3]
        #CIV_fit,object_dictionary['CIV_complex']=line_fitter(sp, "CIV", i,guesses_CIV, limits_CIV, limited_CIV, tied_CIV, xminc, xmaxc,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['CIV_complex'],exclude=exclude, excluding=use_exclude_files,  do_fit=fit_CIV)
        
        CIV_fit1,object_dictionary1['CIV_complex']=line_fitter(sp, "CIV", i,guesses_CIV, limits_CIV, limited_CIV, tied_CIV, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['CIV_complex'], excluding=use_exclude_files,  do_fit=fit_CIV)
        
        

        chi2_1=sp.specfit.chi2/sp.specfit.dof
        dof_1=sp.specfit.dof
        f=copy.data/continuous
        dAIt=(1.0-f)
        dAIc=(1.0-f)[arglim:arglim0]
        dAIw=(1.0-f)[arglim1:arglim2]
        dAIww=(1.0-f)[arglim7:arglim8]
        
        threshold1=np.percentile(dAIt,97.0)
        
        wAIc=dAIc>threshold1
        wAIw=dAIw>threshold1
        wAIww=dAIww>threshold1
        
        
        nabs=1.0*len(dAIc[wAIc])
        nblue=1.0*len(dAIc)

        try:
            rat_flux=nabs/nblue # fraction of points in the lower percentile of the relative  flux (normalized wr to the continuum) between 0 and -10000km/s w/r to the line peak
        except:
            rat_flux=0
        nabs=1.0*len(dAIw[wAIw])
        nblue=1.0*len(dAIw)

        try:
            rat_fluxw=nabs/nblue # fraction of points in the lower percentile of the relative  flux (normalized wr to the continuum) between -2000 and -20000km/s w/r to the line peak
        except:
            rat_fluxw=0
       
        nabs=1.0*len(dAIww[wAIww])
        nblue=1.0*len(dAIww)

        try:
            rat_fluxww=nabs/nblue # fraction of points in the lower percentile of the relative  flux (normalized wr to the continuum) between -2000 and -20000km/s w/r to the line peak
        except:
            rat_fluxww=0
       
     

        print 'flux fraction core=',rat_flux,'nobject=',i
        print  'fux fraction wings=',rat_fluxw
        #BI=
        chi2=chi2_1
        dof=dof_1
        CIV_fit=CIV_fit1
        object_dictionary['CIV_complex']=object_dictionary1['CIV_complex']
        
        
        print 'chi2=',chi2_1
        
        #if chi2>2.0:
        #    continue
        sp=copy.copy()
        sp1=copy.copy()
        
        for j in range(1,4):
            #if use_exclude_files:            
            #    continue
            sp1=copy.copy()
            residuals=sp1.data-(continuous+CIV_fit)
            
            
            delta=treshold
            d=(1.0-delta)**(j)
            

            threshold=np.percentile(residuals,(1-d)*100) 
            
            wlow=residuals<threshold #finding the most negative residuals. 3 iterations each time removing (1-0.97^j)*100 percent of the points 
            sp1.error[wlow]=np.inf
             
            #-----------------testing lines-------------------------#
            #backup=sp1.copy()
            #backup.plotter(xmin=xminc1,xmax=xmaxc1)
            #backup.baseline.powerlaw=False
            #backup.baseline(xmin=xminc1, xmax=xmaxc1, exclude=exclude_cont, subtract=False, reset_selection=False, highlight_fitregion=False,powerlaw=False,quiet=True,LoudDebug=False,annotate=False)
            #continuous=backup.baseline.basespec
            continuous_FUV=continuous
            
            #-----------------testing lines-------------------------#

            
            sp1.data=sp1.data - continuous
            
            
            #sp1 = pyspeckit.Spectrum(xarr=sp.xarr[wlow],data=sp.data[wlow],error=sp.error[wlow])
            

            
            CIV_fit1,object_dictionary1['CIV_complex']=line_fitter(sp1, "CIV", i,guesses_CIV, limits_CIV, limited_CIV, tied_CIV, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['CIV_complex'], excluding=use_exclude_files,  do_fit=fit_CIV)
            chi2=sp1.specfit.chi2/sp1.specfit.dof
            dof=sp1.specfit.dof

            
            if chi2_1<chi2:
                chi2=chi2_1
                dof=dof_1
                break
            chi2_1=chi2
            dof_1=dof
            CIV_fit=CIV_fit1
            object_dictionary['CIV_complex']=object_dictionary1['CIV_complex']
            print 'chi2=',chi2
        sp=copy.copy()
        residuals=sp.data-(continuous+CIV_fit)
        resibal=residuals[arglim:arglim0]
        resicore=residuals[arglim3:arglim4]
        #resicore1=residuals[arglim5:arglim6]
        ratios=sp.data/(continuous+CIV_fit)
        ratcore=ratios[arglim5:arglim6]
        
        threshold1=np.percentile(residuals,3.0)
        threshold2=np.percentile(residuals,97.0)
        threshold3=np.percentile(ratios,3.0)
        
        wlowbal=resibal<threshold1
        wlowcore=ratcore<threshold3
        wupcore=resicore>threshold2
        nabs=1.0*len(resibal[wlowbal])
        nabscore=1.0*len(ratcore[wlowcore])
        nupcore=1.0*len(resicore[wupcore])
        ncore=1.0*len(resicore)
        nblue=1.0*len(resibal)
        ncore1=1.0*len(ratios)
        
        try:
            rat_absw=nabs/nblue # fraction of points in the lowest 3 percentile of the residuals (most negative residualts) in the blue wing between 0 and -10000km/s w/r to the line peak
        except:
            rat_absw=0
        try:
            rat_res=nupcore/ncore # fraction of points in the higheest 3 percentile of the residuals(most positive residuals) in the line core between -1000km and -1000km/s w/r to the line peak
        except:
            rat_res=0
        try:
            rat_absc=nabscore/ncore1
        except:
            rat_absc=0
        print 'abs fraction=',rat_absw,'nobject=',i
        print 'res fraction=',rat_res    
        print 'abs core fraction=',rat_absc,'nobject=',i
        
        #if use_exclude_files==0:            
            
        
        
        sp1=copy.copy()
        sp1.error[wlow]=np.inf
        #-----------------testing lines-------------------------#
        backup=sp1.copy()
        #backup.plotter(xmin=xminc1,xmax=xmaxc1)
        #backup.baseline.powerlaw=False
        #backup.baseline(xmin=xminc1, xmax=xmaxc1, exclude=exclude_cont, subtract=False, reset_selection=False, highlight_fitregion=False,powerlaw=False,quiet=True,LoudDebug=False,annotate=False)
        #continuous=backup.baseline.basespec
        continuous_FUV=continuous
            
        #-----------------testing lines-------------------------#
        

        sp1.data=sp.data - continuous
        
        
        
        CIV_fit1,object_dictionary1['CIV_complex']=line_fitter(sp1, "CIV", i,guesses_CIV, limits_CIV, limited_CIV, tied_CIV, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['CIV_complex'], excluding=use_exclude_files,  do_fit=fit_CIV)
        chi2_1=sp1.specfit.chi2/sp1.specfit.dof
        dof_1=sp1.specfit.dof
        
        if chi2_1<chi2:
            chi2=chi2_1
            dof=dof_1
            CIV_fit=CIV_fit1
            object_dictionary['CIV_complex']=object_dictionary1['CIV_complex']
            
    if line_to_fit=='CIII':
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
        
        CIII_fit,object_dictionary['CIII_complex']=line_fitter(sp, "CIII", i,guesses_CIII, limits_CIII, limited_CIII, tied_CIII, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['CIII_complex'], excluding=use_exclude_files,  do_fit=fit_CIII)
        
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
        print 'chi2=',chi2
        print 'nobject=',i
        
    if line_to_fit=='SiOIV':
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
        
        
        
        SiOIV_fit1,object_dictionary1['SiOIV_complex']=line_fitter(sp, "SiOIV1", i,guesses_SiOIV1, limits_SiOIV1, limited_SiOIV1, tied_SiOIV1, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['SiOIV_complex'], excluding=use_exclude_files,  do_fit=fit_Si)
        
        

        chi2_1=sp.specfit.chi2/sp.specfit.dof
        dof_1=sp.specfit.dof
        f=copy.data/continuous
        dAIt=(1.0-f)
        dAIc=(1.0-f)[arglim:arglim0]
        dAIw=(1.0-f)[arglim1:arglim2]
        threshold1=np.percentile(dAIt,97.0)
        
        wAIc=dAIc>threshold1
        wAIw=dAIw>threshold1
        
        
        nabs=1.0*len(dAIc[wAIc])
        nblue=1.0*len(dAIc)

        try:
            rat_flux=nabs/nblue # fraction of points in the lower percentile of the relative  flux (normalized wr to the continuum) between 0 and -10000km/s w/r to the line peak
        except:
            rat_flux=0
        nabs=1.0*len(dAIw[wAIw])
        nblue=1.0*len(dAIw)

        try:
            rat_fluxw=nabs/nblue # fraction of points in the lower percentile of the relative  flux (normalized wr to the continuum) between -2000 and -20000km/s w/r to the line peak
        except:
            rat_fluxw=0


        print 'flux fraction core=',rat_flux,'nobject=',i
        print  'fux fraction wings=',rat_fluxw
        #BI=
        chi2=chi2_1
        dof=dof_1
        SiOIV_fit=SiOIV_fit1
        object_dictionary['SiOIV_complex']=object_dictionary1['SiOIV_complex']
        
        
        print 'chi2=',chi2_1
        
        #if chi2>2.0:
        #    continue
        sp=copy.copy()
        sp1=copy.copy()
        

        for j in range(1,4):
            
            sp1=copy.copy()
            residuals=sp1.data-(continuous+SiOIV_fit)
            
            
            delta=0.03
            d=(1.0-delta)**(j)
            

            threshold=np.percentile(residuals,(1-d)*100) 
            
            wlow=residuals<threshold #finding the most negative residuals. 3 iterations each time removing (1-0.97^j)*100 percent of the points 
            
            sp1.data=sp1.data - continuous
            sp1.error[wlow]=np.inf
            
            #sp1 = pyspeckit.Spectrum(xarr=sp.xarr[wlow],data=sp.data[wlow],error=sp.error[wlow])
            SiOIV_fit1,object_dictionary1['SiOIV_complex']=line_fitter(sp1,"SiOIV1", i,guesses_SiOIV1, limits_SiOIV1, limited_SiOIV1, tied_SiOIV1, xminc1, xmaxc1, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['SiOIV_complex'], do_fit=fit_Si)

            chi2=sp1.specfit.chi2/sp1.specfit.dof
            dof=sp1.specfit.dof

            
            if chi2_1<chi2:
                chi2=chi2_1
                dof=dof_1
                break
            chi2_1=chi2
            dof_1=dof
            SiOIV_fit=SiOIV_fit1
            
            object_dictionary['SiOIV_complex']=object_dictionary1['SiOIV_complex']
            print 'chi2=',chi2
        sp=copy.copy()
        residuals=sp.data-(continuous+SiOIV_fit)
        resibal=residuals[arglim:arglim0]
        resicore=residuals[arglim3:arglim4]
        
        threshold1=np.percentile(residuals,3.0)
        threshold2=np.percentile(residuals,97.0)
        wlowbal=resibal<threshold1
        wupcore=resicore>threshold2
        nabs=1.0*len(resibal[wlowbal])
        nupcore=1.0*len(resicore[wupcore])
        ncore=1.0*len(resicore)
        nblue=1.0*len(resibal)
        
        
        try:
            rat_absw=nabs/nblue # fraction of points in the lowest 3 percentile of the residuals (most negative residualts) in the blue wing between 0 and -10000km/s w/r to the line peak
        except:
            rat_absw=0
        try:
            rat_res=nupcore/ncore # fraction of points in the higheest 3 percentile of the residuals(most positive residuals) in the line core between -1000km and -1000km/s w/r to the line peak
        except:
            rat_res=0
        print 'abs fraction=',rat_absw,'nobject=',i
        print 'res fraction=',rat_res    
        
        
        sp1.data=sp.data - continuous
        sp1.error[wlow]=np.inf
        
        

        SiOIV_fit,object_dictionary['SiOIV_complex']=line_fitter(sp1,"SiOIV1", i,guesses_SiOIV1, limits_SiOIV1, limited_SiOIV1, tied_SiOIV1, xminc1, xmaxc1, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['SiOIV_complex'], do_fit=fit_Si)
        chi2_1=sp1.specfit.chi2/sp1.specfit.dof
        dof_1=sp1.specfit.dof
        
        if chi2_1<chi2:
            chi2=chi2_1
            dof=dof_1
            SiOIV_fit=SiOIV_fit1
            object_dictionary['SiOIV_complex']=object_dictionary1['SiOIV_complex']
            
















        """

        SiOIV_fit,object_dictionary['SiOIV_complex']=line_fitter(sp,"SiOIV1", i,guesses_SiOIV1, limits_SiOIV1, limited_SiOIV1, tied_SiOIV1, xminc, xmaxc, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['SiOIV_complex'], do_fit=fit_Si)
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
        
        print 'chi2=',chi2
        if chi2>2.0:
            continue
        """
    if line_to_fit=='Lya':
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
        
        Lya_fit,object_dictionary['Lya_complex']=line_fitter(sp,"Lya", i,guesses_Lya, limits_Lya, limited_Lya, tied_Lya, xminc1, xmaxc1, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['Lya_complex'], do_fit=fit_Si)
        
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
        print 'chi2=',chi2
        
    if line_to_fit=='C':        
        
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
        
        C_fit,object_dictionary['C_complex']=line_fitter(sp, "C", i,guesses_C, limits_C, limited_C, tied_C, xminc1, xmaxc1,magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['C_complex'], excluding=use_exclude_files,  do_fit=fit_C)
        
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
        print 'chi2=',chi2
    if line_to_fit=='Halpha':
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]


        Halpha_fit,object_dictionary['Halpha_complex']=line_fitter(sp,"Halpha", i, guesses_Halpha, limits_Halpha, limited_Halpha, tied_Halpha,
                                                                   xmin_Halpha, xmax_Halpha, magorder=mag_order,plot_path=plot_objpath,
                                                                   linenames=lines_dict['Halpha_complex'],  do_fit=fit_Halpha)
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
        print 'chi2=',chi2

    
        
    #del(number,mag,redshift)
    
    if fileroot in sample_dictionary:
        try:
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
            print 'opened'
        except:
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']={}
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']={}
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
        
        print fileroot, i
    else:
        print fileroot, ' does not exists'
        sample_dictionary[fileroot]={}
        sample_dictionary[fileroot]['model']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
        
        

    

    
    
        
    sp = pyspeckit.Spectrum(spectrum_file)
    sp.crop(argxmin,argxmax)
    #sp.data=sp.data[~wlow]
    #sp.error=sp.error[~wlow]
    #sp.xarr=sp.xarr[~wlow]
    #continuous=continuous[~wlow]
    # -----------set up unit properly------------------#        
    mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    #sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'

    sp.data = 10**(-1.0*mag_order)*np.array(sp.data.tolist())    
    sp.error= 10**(-1.0*mag_order)*np.array(sp.error.tolist())    
    #-------------- set up unit properly------------#
    
    #----Superposition plot --------#
    total=(continuous +  CIV_fit + C_fit +
           + SiOIV_fit +CIII_fit+Lya_fit+Halpha_fit)
    total_lines=( CIV_fit + C_fit +
           + SiOIV_fit +CIII_fit+Lya_fit+Halpha_fit)
    
    plot_file=plot_objpath + line_to_fit + "_"+ fileroot.split('.')[0] + ".png"
    
    if plot_best_fits:
    
        pylab.rcParams["figure.figsize"]=16,6
        copy1=sp.copy() 

        argxmin=np.argmin(np.abs(sp.xarr.value-xminc))
        argxmax=np.argmin(np.abs(sp.xarr.value-xmaxc))
        copy1.crop(argxmin,argxmax)
    
        pylab.figure()
        pylab.ylim(ymin=0,ymax=1.1*copy1.data.max())
        
        pylab.xlim(xmin=xminc-100.0,xmax=xmaxc+100.0)
        #pylab.ylabel(r'$10^{'+str(mag_order-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        #pylab.xlabel(r'$\AA$')
        pylab.yscale('linear')
        try:
            pylab.plot(sp.xarr,sp.data,'k',label=' fluxcenter rat='+str(np.round(rat_flux,2))
                       +' fluxwingfar rat='+str(np.round(rat_fluxww,2)) +' fluxwing rat='+str(np.round(rat_fluxw,2)) + '  SN='+str(np.round(SN,2))+'chi2='+str(np.round(chi2,2))+'  abs rat='+str(np.round(rat_absw,2))+'  exccess rat='+str(np.round(rat_res,2)) +'  abs rat core='+str(np.round(rat_absc,2)))
        
        except:
            pylab.plot(sp.xarr,sp.data,'k')

        pylab.plot(sp.xarr, total,'r')
        pylab.plot(sp.xarr, total_lines,'r')
        total=continuous
        pylab.plot(sp.xarr, total,'gray')
        #pylab.legend(fontsize=10)
        #sp.plotter.figure.savefig(plot_file,format='pdf', dpi=600,bbox_inches='tight')
        pylab.savefig(plot_file)
    
        pylab.close('all')
        #----Superposition plot --------#



    



    
    obj_path=dataout_path+fileroot + '/'
    
    if not os.path.exists(obj_path):
        os.mkdir(obj_path)
    np.savetxt(obj_path +'obs_spectrum.txt' ,np.transpose([sp.xarr,sp.data]))
        
    
    if line_to_fit=='CIV':
        np.savetxt(obj_path +'CIV_fit.txt'          ,np.transpose([sp.xarr,CIV_fit]) )
        np.savetxt(obj_path +'continuous_CIV.txt' ,np.transpose([sp.xarr,continuous]))
        
        
                
    if line_to_fit=='CIII':
        np.savetxt(obj_path +'CIII_fit.txt'          ,np.transpose([sp.xarr,CIII_fit]) )
        np.savetxt(obj_path +'continuous_CIII.txt' ,np.transpose([sp.xarr,continuous]))
        
    
    if line_to_fit=='SiOIV':
        np.savetxt(obj_path +'SiOIV_fit.txt'         ,np.transpose([sp.xarr,SiOIV_fit]) )
        np.savetxt(obj_path +'continuous_SiOIV.txt' ,np.transpose([sp.xarr,continuous]))
    
    if line_to_fit=='Lya':
        np.savetxt(obj_path +'Lya_fit.txt'         ,np.transpose([sp.xarr,Lya_fit]) )
        np.savetxt(obj_path +'continuous_Lya.txt' ,np.transpose([sp.xarr,continuous]))
    
    if line_to_fit=='Halpha':
        np.savetxt(obj_path +'Halpha_fit.txt'         ,np.transpose([sp.xarr,Halpha_fit]) )
        np.savetxt(obj_path +'continuous_Halpha.txt' ,np.transpose([sp.xarr,continuous]))
    
        
    
    
    
    #np.savetxt(obj_path +'H_fit.txt'     ,np.transpose([sp.xarr,H_fit]) )
    
    sample_dictionary[fileroot]['model']    
    


    
    
    

    if line_to_fit=='CIV':
        sample_dictionary[fileroot]['model']['CIV_complex']['datafile']=         obj_path +'CIV_fit.txt'
        sample_dictionary[fileroot]['model']['CIV_complex']['continuous']={}
        sample_dictionary[fileroot]['model']['CIV_complex']['continuous']['datafile']=        obj_path + 'continuous_CIV.txt'
        sample_dictionary[fileroot]['model']['CIV_complex']['chi2']=        chi2
        sample_dictionary[fileroot]['model']['CIV_complex']['dof']=        dof
        
        sample_dictionary[fileroot]['model']['CIV_complex']['abs_ratio']=rat_absw
        sample_dictionary[fileroot]['model']['CIV_complex']['abs_core_ratio']=rat_absc
        
        
        
        sample_dictionary[fileroot]['model']['CIV_complex']['res_ratio']=rat_res
        sample_dictionary[fileroot]['model']['CIV_complex']['core_ratio']=rat_flux 
        sample_dictionary[fileroot]['model']['CIV_complex']['wing_ratio']=rat_fluxw
        sample_dictionary[fileroot]['model']['CIV_complex']['wing_far_ratio']=rat_fluxww
        sample_dictionary[fileroot]['model']['CIV_complex']['SN']=        SN
        
    if line_to_fit=='CIII':
        sample_dictionary[fileroot]['model']['CIII_complex']['datafile']=         obj_path +'CIII_fit.txt'
        sample_dictionary[fileroot]['model']['CIII_complex']['continuous']={}
        sample_dictionary[fileroot]['model']['CIII_complex']['continuous']['datafile']=        obj_path + 'continuous_CIII.txt'
        sample_dictionary[fileroot]['model']['CIII_complex']['continuous']['chi2']=        chi2
        sample_dictionary[fileroot]['model']['CIII_complex']['continuous']['dof']=        dof
    if line_to_fit=='SiOIV':           
        sample_dictionary[fileroot]['model']['SiOIV_complex']['datafile']=      obj_path +'SiOIV_fit.txt'
        sample_dictionary[fileroot]['model']['SiOIV_complex']['continuous']={}
        sample_dictionary[fileroot]['model']['SiOIV_complex']['continuous']['datafile']=        obj_path + 'continuous_SiOIV.txt'
        sample_dictionary[fileroot]['model']['SiOIV_complex']['continuous']['chi2']=        chi2
        sample_dictionary[fileroot]['model']['SiOIV_complex']['continuous']['dof']=        dof

        
        sample_dictionary[fileroot]['model']['SiOIV_complex']['abs_ratio']=rat_absw
        sample_dictionary[fileroot]['model']['SiOIV_complex']['res_ratio']=rat_res
        sample_dictionary[fileroot]['model']['SiOIV_complex']['core_ratio']=rat_flux 
        sample_dictionary[fileroot]['model']['SiOIV_complex']['wing_ratio']=rat_fluxw
        sample_dictionary[fileroot]['model']['SiOIV_complex']['SN']=        SN
    if line_to_fit=='Lya':           
        sample_dictionary[fileroot]['model']['Lya_complex']['datafile']=      obj_path +'Lya_fit.txt'
        sample_dictionary[fileroot]['model']['Lya_complex']['continuous']={}
        sample_dictionary[fileroot]['model']['Lya_complex']['continuous']['datafile']=        obj_path + 'continuous_Lya.txt'
        sample_dictionary[fileroot]['model']['Lya_complex']['continuous']['chi2']=        chi2
        sample_dictionary[fileroot]['model']['Lya_complex']['continuous']['dof']=        dof 
    if line_to_fit=='Halpha':           
        sample_dictionary[fileroot]['model']['Halpha_complex']['datafile']=      obj_path +'Lya_fit.txt'
        sample_dictionary[fileroot]['model']['Halpha_complex']['continuous']={}
        sample_dictionary[fileroot]['model']['Halpha_complex']['continuous']['datafile']=        obj_path + 'continuous_Halpha.txt'
        sample_dictionary[fileroot]['model']['Halpha_complex']['continuous']['chi2']=        chi2
        sample_dictionary[fileroot]['model']['Halpha_complex']['continuous']['dof']=        dof
       
        

    sample_dictionary[fileroot]['obs_spectrum']= obj_path +'obs_spectrum.txt'
    
    with open(obj_path + line_to_fit +'_dictionary.json', 'wb') as fp:
        #json.dump(object_dictionary[line_to_fit+'_complex'], fp)
        json.dump(sample_dictionary[fileroot]['model'][line_to_fit+'_complex'], fp)
    with open(dataout_path + 'sample_dictionary.json', 'wb') as fp:
        json.dump(sample_dictionary, fp)


    t1=time.time()
    print 'elapsed time=', t1-t0
with open(dataout_path + 'sample_dictionary.json', 'wb') as fp:
    json.dump(sample_dictionary, fp)


#with open('data.json', 'rb') as fp:
#    data = json.load(fp)
