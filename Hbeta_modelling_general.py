# to run python Hbeta_modelling_general.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name # remember to after 'path_to_data' to include the symbol '/'
import matplotlib
#matplotlib.use('Agg') 
matplotlib.use('Agg')
import time
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


plot_best_fits=1 # 1 if you want to plot the best fits. to run in cluster must be set to 0
treshold=0.05 #this is the percentage of most negative pixels that are going to be exclude in the fit after the continuum removal
              # this procedure is done iteratively for 3 times.
              #For low-z quasars with little absortion treshold=0.03 is ok (z<3). For  high-z quasars treshold=0.05 would be better
use_exclude_files=0 # 0 if you want the code to fit the lines witouth using manually selected continuum windwos and manually selected 
                    #absorption features. The files should be in the folder ./excludes/ with the following filenames:
                    # exclude_cont_CIV.txt and exclude_CIV.txt. Replace CIV for the name of other lines
                    # you want to fit e.g. CIII, SiOIV, Lya. 

extension='txt' # please provide the common extension of the spectrum filenames

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
    line_to_fit='Hbeta'
else:
    line_to_fit=sys.argv[1]
    
    

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
    





line_to_fit='Hbeta'

rcParams['text.usetex'] = True #Very important to force python to recognize Latex code
Hbeta_iter=3
final_Hb_fit=1
#------Executing the configuration file that includes the requiered parameters---#
execfile("./constraints_single.cfg") 




xmin=4400
xmax=5600
n_spec=1










xminc1=4400.0
xmaxc1=5300.0
#exclude_cont=[100,2640,2670,3020,3040,10000]

exclude_cont=[100.0,   4400.0,    4450,    5400.2,    5450.1,    10500]
xminc1=exclude_cont[1]
xmaxc1=exclude_cont[-2]




if use_exclude_files:
    exclude_file = "./excludes/exclude_cont_Hbeta.txt"
    exclude_conti=np.loadtxt(exclude_file,skiprows=2)
    exclude_cont=[100,2640,2670,3020,3040,10000]
    xminc1=exclude_cont[1]
    xmaxc1=exclude_cont[-2]

                












try:
    sp_to_run=range(len(filenames))
except:
    filenames=np.array([str(filenames)])
    print 'only one file'
    sp_to_run=[0]


#sp_to_run=range(len(plate))[20:100]

#indexes_try=np.append(indexnoabs[:50],indexbal[:50])            
#indexes_try=np.sort(indexes_try)
#sp_to_run=indexes_try[:10]


for i in sp_to_run:

    
    t0=time.time()
    
    fileroot= filenames[i]
    objname=fileroot.split('.'+extension)[0]
    
    spectrum_file=file=data_dir+fileroot
    print 'fitting ',line_to_fit,' in ', fileroot,'\n\n\n\n'


    object_dictionary={}
    object_dictionary1={}
    

        
    plot_objpath=plots_path
    if not os.path.exists(plot_objpath):
        os.mkdir(plot_objpath)

    
    

    
    sp = pyspeckit.Spectrum(spectrum_file)
    
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
    
    fe_template_OP=sp.copy()
    fe_template_OP.data=0.0*sp.data
    continuous=0.0*sp.data 
    Halpha_fit=0.0*sp.data 
    Hbeta_fit=0.0*sp.data 
    C_fit=0.0*sp.data 
    fe_fit=0.0*sp.data 
    fe_Hb=0.0*sp.data 
    balmer_data=0.0*sp.data 
    MgII_fit=0.0*sp.data 
    SiOIV_fit=0.0*sp.data 
    OI_fit=0.0*sp.data 
    OII_fit=0.0*sp.data 
    Lya_fit=0.0*sp.data 
    CIISF_fit=0.0*sp.data 
    NIVF_fit=0.0*sp.data 
    AlIIF_fit=0.0*sp.data 
    Hgamma_fit=0.0*sp.data 
    Hdelta_fit=0.0*sp.data 
    balmer_tot=0.0*sp.data 
    

    #------------opening spectrum i---------------------#
    
        
    # -----------set up units properly------------------#        
    try:
        mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    except:
        continue
        #w1=~np.isnan(np.log10(np.mean(sp.data)))
        #w2=~np.isinf(np.log10(np.mean(sp.data)))
        #mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data[w1*w2]))))
    sp.xarr.unit==un.Angstrom

    #sp.xarr.units='angstrom'
    sp.xarr.xtype = 'wavelength'
    sp.xarr.unit==un.Angstrom
    #sp.xarr.units='Angstrom'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.data = 10**(-1.0*mag_order)*np.array(sp.data.tolist())    
    sp.error= 10**(-1.0*mag_order)*np.array(sp.error.tolist())    
    #-------------- set up unit properly------------#
    copy=sp.copy()
        
    #-------------continuous fitting-------------------#
    wlmin=sp.xarr.value[0]
    wlmax=sp.xarr.value[-1]

    
    









    
    sp1=sp.copy()
    fe_Hb=0.0*sp.data
    f_name=["rescaling0","rescaling1","rescaling2","rescaling3","rescaling4"]
    if line_to_fit=='Hbeta':
        for iter in range(Hbeta_iter):
            if use_exclude_files:
                exclude_cont=np.loadtxt(exclude_file,skiprows=2)
                exclude_cont=exclude_cont[i,:][:]
                xminc1=exclude_cont[1]
                xmaxc1=exclude_cont[-2]
            
            fe_template,galaxy_template=load_op_template(sp,mag_order,fe_template_OP_file,galaxy_template_file)
            
            
            backup=sp.copy()
            #backup.plotter()
            backup.baseline.powerlaw=False
            backup.baseline(xmin=xminc1, xmax=xmaxc1, exclude=exclude_cont, subtract=False, reset_selection=False, highlight_fitregion=False,powerlaw=False,interactive=False,qxuiet=False,LoudDebug=False,annotate=False)
            continuous=backup.baseline.basespec
            continuous_OP=continuous
            
            cont_OP2=continuous_OP
            xmin_fe=4200
            xmax_fe=5600
            
            
            
            sp_0=sp.data
            sp.data=sp1.data - continuous_OP
            cont=sp.copy()
            cont.data=continuous
            sp_1=sp.data
            if iter<=1:
                Hbeta_fit,object_dictionary['Hbeta_complex']=line_fitter(sp, "Hbeta", i, guesses_Hbeta0, limits_Hbeta0, limited_Hbeta0, tied_Hbeta0, xmin,
                                                                         xmax,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['Hbeta_complex'],  do_fit=fit_Hbeta)
                
            else:
                Hbeta_fit,object_dictionary['Hbeta_complex']=line_fitter(sp, "Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin,
                                                                         xmax,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['Hbeta_complex'],  do_fit=fit_Hbeta)
                
                
            fwhm1=wl_to_kms(object_dictionary['Hbeta_complex']['Hbeta1']['fwhm'],4682)
            fwhm2=wl_to_kms(object_dictionary['Hbeta_complex']['Hbeta2']['fwhm'],4682)
            fwhm_low=np.amin([fwhm1,fwhm2])
            if iter>0:
                fwhm_low=np.amin([fwhm_low,fwhm_low1])
            fwhm_low1=fwhm_low
                
            #pylab.ylim(-1.0,1.0)
            
            
            print fwhm1, fwhm2, fwhm_low
            
            fe_op_params=np.ones(3)
            fe_op_params[0]=fwhm_low
            fe_op_params[1]=0.0
            
            lambda0=4860.0
            fe_Hb,cont,fe_template_OP,fe_params,chi2op=fe_scale(sp,fe_op_params,"optical",i,xmin_fe,xmax_fe,fe_template,cont,galaxy_template,f_name[iter],mag_order,lambda0,objectname=objname,plot_path=plot_objpath, do_fit=fit_fe_OP)
                
                
            
            pylab.rcParams["figure.figsize"]=16,6
            xmin_fe=4200.0
            xmax_fe=5600.0
            print "xmin,xmax=   ",xmin_fe,xmax_fe
                
                
            sp.data=sp1.data-fe_Hb
            #if(iter==2):
            #    sys.exit(1)

            print fe_op_params
        print "fe_op_params=",fe_op_params
            
        fe_op_params[2]=fe_params[0]
        fe_op_params[1]+=fe_params[1]
        print "fe_op_params1=",fe_op_params
        
        continuous_OP=cont.data
        sp.data=sp.data-continuous_OP
        print np.sum(galaxy_template.data)
        lambda0=4860.0
        #fe_Hb,cont,fe_template_OP,fe_params=fe_scale(sp,fe_op_params,"optical",i,xmin_fe,xmax_fe,fe_template,cont,galaxy_template,mag_order,lambda0,plot_path=plot_objpath,do_fit=fit_fe_OP)
        print np.sum(galaxy_template.data)
            
            
            

        continuous_OP=cont.data #+ galaxy_template.data
        
        continuous=continuous_OP
        pylab.rcParams["figure.figsize"]=8,6
        
        Hgamma_fit=0.0*sp.data
        Hdelta_fit=0.0*sp.data
        if final_Hb_fit==1:
            Hbeta_fit,object_dictionary['Hbeta_complex']=line_fitter(sp, "Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin_Hbeta,
                                                                             xmax_Hbeta,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['Hbeta_complex'],  do_fit=fit_Hbeta)




        sp.data=sp1.data
        argfe_max=np.argmin(np.abs(sp.xarr-fe_template_OP.xarr[-1]))
        argfe_min=np.argmin(np.abs(sp.xarr-fe_template_OP.xarr[0]))
        if sp.xarr[argfe_max]<fe_template_OP.xarr[-1]: argfe_max+=1
        if sp.xarr[argfe_min]>fe_template_OP.xarr[-1]: argfe_min-=1
        new_fex=np.concatenate([ sp.xarr[:argfe_min+1], fe_template_OP.xarr , sp.xarr[argfe_max:]   ] )
        new_fey=np.zeros_like(new_fex)
        new_fey[argfe_min:argfe_min+len(fe_template_OP.xarr)]=fe_template_OP.data
                #fe_Hb=new_fey
            
            
        continuous1=np.interp(new_fex,sp.xarr,continuous)
        fe_fit1=np.interp(new_fex,sp.xarr,fe_fit)
        balmer_data1=np.interp(new_fex,sp.xarr,balmer_data)

    

        print 'done with nobject=',i, fileroot
        chi2=sp.specfit.chi2/sp.specfit.dof
        dof=sp.specfit.dof
    
        
    #del(number,mag,redshift)
    
    if fileroot in sample_dictionary:
        try:
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
        except:
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']={}
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']={}
            sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
        
        print fileroot, i
    else:
        sample_dictionary[fileroot]={}
        sample_dictionary[fileroot]['model']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']={}
        sample_dictionary[fileroot]['model'][line_to_fit+'_complex']['lines']=object_dictionary[line_to_fit+'_complex']
        
        


    

    
    
        
    sp = pyspeckit.Spectrum(spectrum_file)
    try:
        sp.crop(argxmin,argxmax)
    except:
        continue
    
    
    #sp.data=sp.data[~wlow]
    #sp.error=sp.error[~wlow]
    #sp.xarr=sp.xarr[~wlow]
    #continuous=continuous[~wlow]
    # -----------set up unit properly------------------#        
    mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    #sp.xarr.units='Angstroms'
    sp.xarr.unit==un.Angstrom
    sp.xarr.xtype = 'wavelength'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    #sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'

    sp.data = 10**(-1.0*mag_order)*np.array(sp.data.tolist())    
    sp.error= 10**(-1.0*mag_order)*np.array(sp.error.tolist())    
    #-------------- set up unit properly------------#
    
    #----Superposition plot --------#
    

    total=(continuous +  Halpha_fit + Hbeta_fit + C_fit +
           fe_fit  + fe_Hb + balmer_data +  MgII_fit 
           + SiOIV_fit +CIISF_fit+ NIVF_fit + AlIIF_fit+
           Hgamma_fit + Hdelta_fit+OI_fit+OII_fit)
    total_lines=(Halpha_fit + Hbeta_fit + C_fit +
                 MgII_fit 
                 + SiOIV_fit +CIISF_fit+ NIVF_fit + AlIIF_fit+
                 Hgamma_fit + Hdelta_fit+OI_fit+OII_fit)


    if plot_best_fits:
        plot_file=plot_objpath + line_to_fit + "_"+ fileroot.split('.')[0] + ".png"

    
        pylab.rcParams["figure.figsize"]=16,6
        copy1=sp.copy() 
        
        argxmin=np.argmin(np.abs(sp.xarr.value-xmin))
        argxmax=np.argmin(np.abs(sp.xarr.value-xmax))
        copy1.crop(argxmin,argxmax)
    
        pylab.figure()
        pylab.ylim(ymin=0,ymax=1.1*copy1.data.max())
        
        pylab.xlim(xmin=xmin,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.yscale('linear')
        try:
            pylab.plot(sp.xarr.value,sp.data,'k',label=' fluxcenter rat='+str(np.round(rat_flux,2))
                       +' fluxwingfar rat='+str(np.round(rat_fluxww,2)) +' fluxwing rat='+str(np.round(rat_fluxw,2)) + '  SN='+str(np.round(SN,2))+'chi2='+str(np.round(chi2,2))+'  abs rat='+str(np.round(rat_absw,2))+'  exccess rat='+str(np.round(rat_res,2)) +'  abs rat core='+str(np.round(rat_absc,2)))
        
        except:
            pylab.plot(sp.xarr.value,sp.data,'k')

            pylab.plot(sp.xarr.value, total,'r')
            pylab.plot(sp.xarr.value, total_lines,'r')
            
            total=continuous + fe_fit  + balmer_data + fe_Hb
            pylab.plot(sp.xarr, total,'b')
            pylab.plot(sp.xarr, fe_fit+fe_Hb,'b')
            

            total=continuous
            pylab.plot(sp.xarr.value, total,'gray')
            pylab.legend(fontsize=10)
            #sp.plotter.figure.savefig(plot_file,format='png', dpi=600,bbox_inches='tight')
            pylab.savefig(plot_file,format='png', dpi=600,bbox_inches='tight')
            pylab.close('all')
    #----Superposition plot --------#



    



    
    obj_path=dataout_path+fileroot + '/'
    
    if not os.path.exists(obj_path):
        os.mkdir(obj_path)
    np.savetxt(obj_path +'obs_spectrum.txt' ,np.transpose([sp.xarr.value,sp.data]))
        
    
    np.savetxt(obj_path +'Hbeta_fit.txt'         ,np.transpose([sp.xarr.value,Hbeta_fit]) )
    np.savetxt(obj_path +'feOp_fit.txt'          ,np.transpose([sp.xarr.value,fe_fit+fe_Hb]) )
    np.savetxt(obj_path +'continuous_Hbeta.txt' ,np.transpose([sp.xarr.value,continuous]))
    
    
    
    #np.savetxt(obj_path +'H_fit.txt'     ,np.transpose([sp.xarr.value,H_fit]) )
    
    sample_dictionary[fileroot]['model']    
    

    sample_dictionary[fileroot]['model']['feOp']={}
    

    sample_dictionary[fileroot]['model']['Hbeta_complex']['datafile']=      obj_path +'Hbeta_fit.txt'
    sample_dictionary[fileroot]['model']['feOp']['datafile']=              obj_path +'feOp_fit.txt'
    sample_dictionary[fileroot]['model']['feOp']['fwhm']=                  fe_op_params[0]
    sample_dictionary[fileroot]['model']['feOp']['shift']=                 fe_op_params[1]#    sample_dictionary[fileroot]['model']['feOp']['scale']=                 fe_uv_params[2]
    sample_dictionary[fileroot]['model']['feOp']['datafile']=              obj_path +'feOp_fit.txt'
    sample_dictionary[fileroot]['model']['Hbeta_complex']['continuous']={}
    sample_dictionary[fileroot]['model']['Hbeta_complex']['continuous']['datafile']=        obj_path + 'continuous_Hbeta.txt'
    sample_dictionary[fileroot]['model']['Hbeta_complex']['chi2']=        chi2
    sample_dictionary[fileroot]['model']['Hbeta_complex']['dof']=        dof

    


    sample_dictionary[fileroot]['obs_spectrum']= obj_path +'obs_spectrum.txt'
    
    with open(obj_path + line_to_fit +'_dictionary.json', 'wb') as fp:
        #json.dump(object_dictionary[line_to_fit+'_complex'], fp)
        json.dump(sample_dictionary[fileroot]['model'][line_to_fit+'_complex'], fp)
    with open(dataout_path + 'sample_dictionary.json', 'wb') as fp:
        json.dump(sample_dictionary, fp)


    t1=time.time()
    print 'elapsed time=',t1-t0
with open(dataout_path + 'sample_dictionary.json', 'wb') as fp:
    json.dump(sample_dictionary, fp)


#with open('data.json', 'rb') as fp:
#    data = json.load(fp)
