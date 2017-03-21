# to run python MgII_modelling_general.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name # remember to after 'path_to_data' to include the symbol '/'
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
    line_to_fit='MgII'
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
    





line_to_fit='MgII'

rcParams['text.usetex'] = True #Very important to force python to recognize Latex code
Hbeta_iter=3
final_Hb_fit=1
#------Executing the configuration file that includes the requiered parameters---#
execfile("./constraints_single.cfg") 




xmin=1100
xmax=9000
n_spec=1









xminc =2655.0
xmaxc=3020.0
exclude_cont=[100,2640,2670,3020,3040,10000]
        
xminc1=exclude_cont[1]
xmaxc1=exclude_cont[-2]

if use_exclude_files:
    exclude_file = "./excludes/exclude_cont_MgII.txt"
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
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.unit = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $\AA^{-1}$'
    sp.data = 10**(-1.0*mag_order)*np.array(sp.data.tolist())    
    sp.error= 10**(-1.0*mag_order)*np.array(sp.error.tolist())    
    #-------------- set up unit properly------------#
    copy=sp.copy()
        
    #-------------continuous fitting-------------------#
    wlmin=sp.xarr.value[0]
    wlmax=sp.xarr.value[-1]

    
    wlimin_FUV=xminc


    
            
    sp2=sp.copy()
    fe_fit=0.0*sp.data
    f_name=["rescalinguv0","rescalinguv1","rescalinguv2","rescalinguv3"]
    print 'fitting nobject=',i, fileroot
    for inter in range(3):
        #sp,continuous, balmer_tot,balmer_template,L_model,scale_balmer,index,wlmin_UV=uv_continuum(i,sp,mag_order,UV_limits,w,balmer_cont_template,balmer_lines_template,plot_objpath) 
        #continuous_UV=continuous
        wlimin_UV=xminc
                
                
                


        xminc=exclude_cont[1]
        xmaxc=exclude_cont[-2]
        
        
        if use_exclude_files:            
            exclude_cont=exclude_conti[i,:][:]
            xminc1=exclude_cont[1]
            xmaxc1=exclude_cont[-2]
            print xminc1, xmaxc1, 'cont limits!!!!'

        backup=sp.copy()
        #backup.plotter(xmin=xmin,xmax=6900)
        backup.baseline.powerlaw=False
        backup.baseline(xmin=xminc, xmax=xmaxc, exclude=exclude_cont, subtract=False, reset_selection=False, highlight_fitregion=False,powerlaw=False,quiet=True,LoudDebug=False,annotate=False)
                        
            
        continuous=backup.baseline.basespec
        continuous_UV=continuous
        
        sp.data=sp.data - continuous_UV  + fe_fit 
        sp1=sp.copy()

        fe_template_UV=load_uv_template(sp,mag_order,fe_template_UV_file)
        #-------------Fitting the Blue Bump------------------------
        xmin_fe=2550
        xmax_fe=3020
        if inter==0:
            fe_fit,fe_uv_params=fe_fitter(sp,"SBB",i,xmin_fe,xmax_fe,fe_template_UV,mag_order,plot_path=plot_objpath, do_fit=fit_fe_UV)
            
            
            pylab.rcParams["figure.figsize"]=8,6            
            MgII_fit,object_dictionary['MgII_complex']=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['MgII_complex'], do_fit=fit_Mg)
            #if fit_Mg: print "FWHM_Mg=", sp.specfit.measure_approximate_fwhm()
            try:
                FWHM_Mg=sp.specfit.measure_approximate_fwhm().value
            except:
                print 'FWHM(MgII) cannot be obtained in '+ fileroot.split('.txt')[0]
                FWHM_Mg=np.inf
                break
            FWHMV_Mg=wl_to_kms(FWHM_Mg,2800)
            # SAME IRON FWHM THAN MG 
            sp.data=sp1.data  # DO NOT FORGET TO INCLUDE THIS LINE
            sp.error=sp1.error
            sigma_template=900/2.355
            sigmav_Mg=FWHMV_Mg/2.355
            sigma=kms_to_wl(np.sqrt(sigmav_Mg*sigmav_Mg-sigma_template*sigma_template),2800)  
        
        #---new removing residuals---#    
        residuals=sp1.data-fe_fit-MgII_fit
        threshold=np.percentile(residuals,3.0*(inter+1))
        wlow=residuals<threshold
        sp.error[wlow]=np.inf
        #---new removing residuals---#    


        
        fe_template_UV=load_uv_template(sp,mag_order,fe_template_UV_file)
        fe_template_UV.data=gauss_conv(fe_template_UV.data, sigma, order=0, output=None, mode='reflect', cval=0.0)
        fe_template_UV.data*=fe_uv_params[2]
        
        pylab.rcParams["figure.figsize"]=16,6            
        lambda0=2800.0
        #fe_fit,cont,fe_template_UV,fe_params=fe_scale(sp,fe_uv_params,"SBB",i,xmin_fe,xmax_fe,fe_template_UV,cont,galaxy_template,f_name[inter],mag_order,lambda0,plot_path=plot_objpath, do_fit=fit_fe_UV)
        fe_fit,params1=rescale(sp,"SBB",i,xmin_fe,xmax_fe,fe_template_UV,mag_order,plot_path=plot_objpath, do_fit=fit_fe_UV)
        fe_uv_params[0]=FWHMV_Mg
        fe_uv_params[1]=params1[1]
        fe_uv_params[2]=params1[0]
        pylab.rcParams["figure.figsize"]=8,6            
        MgII_fit,object_dictionary['MgII_complex']=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['MgII_complex'], do_fit=fit_Mg)
        
        try:
            FWHM_Mg=sp.specfit.measure_approximate_fwhm().value
        except:
            print 'FWHM(MgII) cannot be obtained in '+ fileroot.split('.txt')[0]
            FWHM_Mg=np.inf
            break
            
        #FWHM_Mg=sp.specfit.measure_approximate_fwhm().value
        FWHMV_Mg=wl_to_kms(FWHM_Mg,2800)
                
                
        

                
                
        #TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE
        #TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE
            
        sp.data=sp1.data  # DO NOT FORGET TO INCLUDE THIS LINE
        sp.error=sp1.error

        
        #---new removing residuals---#    
        residuals=sp1.data-fe_fit-MgII_fit
        threshold=np.percentile(residuals,3.0*(inter+1))
        wlow=residuals<threshold
        sp.error[wlow]=np.inf
        #---new removing residuals---#    



        sigma_template=900/2.355
        sigmav_Mg=FWHMV_Mg/2.355
        sigma=kms_to_wl(np.sqrt(sigmav_Mg*sigmav_Mg-sigma_template*sigma_template),2800) 
        fe_template_UV.data=gauss_conv(fe_template_UV.data, sigma, order=0, output=None, mode='reflect', cval=0.0)
        fe_template_UV.data*=fe_uv_params[2]
                
        pylab.rcParams["figure.figsize"]=16,6            
        
        fe_fit,params1=rescale(sp,"SBB",i,xmin_fe,xmax_fe,fe_template_UV,mag_order,plot_path=plot_objpath, do_fit=fit_fe_UV)
        fe_uv_params[0]=FWHMV_Mg
        fe_uv_params[1]=params1[1]
        fe_uv_params[2]=params1[0]
        print "fe_uv_params=",fe_uv_params
        pylab.rcParams["figure.figsize"]=8,6            
        MgII_fit,object_dictionary['MgII_complex']=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['MgII_complex'], do_fit=fit_Mg)
        
        fe_uv_params1=fe_uv_params
        
        try:
            FWHM_Mg=sp.specfit.measure_approximate_fwhm().value
        except:
            print 'FWHM(MgII) cannot be obtained in '+ fileroot.split('.txt')[0]
            FWHM_Mg=np.inf
            break
            
        #FWHM_Mg=sp.specfit.measure_approximate_fwhm().value
        FWHMV_Mg=wl_to_kms(FWHM_Mg,2800)
        
        sigma_template=900/2.355
        sigmav_Mg=FWHMV_Mg/2.355
        sigma=kms_to_wl(np.sqrt(sigmav_Mg*sigmav_Mg-sigma_template*sigma_template),2800) 
        
        sp.data=sp2.data -fe_fit
        sp.error=sp1.error
        #TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE
        #TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE#TWICE
    if FWHM_Mg==np.inf:
        continue
    fe_fit1=fe_fit
    sp.data=sp2.data - continuous_UV
    fe_template_UV.xarr=sp.xarr
    fe_template_UV.data=fe_fit1
    
    
            
            
    xbreak=2810.0
    xbalmer=3648.0
    fitmin=np.argmin(np.abs(sp.xarr.value-xmin_fe))
    fitmax=np.argmin(np.abs(sp.xarr.value-xmax_fe))
            
            
    sp.data=sp2.data - continuous_UV
    #fe_fit_red,params1=rescale(sp,"SBB",i,xbreak,xmax_fe,fe_template_UV,mag_order,plot_path=plot_objpath,slim=[0.92,1.02], vlim=[-0.1,0.1], do_fit=fit_fe_UV)
    fe_fit_red,params1=rescale(sp,"SBB",i,xbreak,xmax_fe,fe_template_UV,mag_order,plot_path=plot_objpath,slim=[0.70,1.10], vlim=[-0.1,0.1], do_fit=fit_fe_UV)
    
    
    fe_template_UV.data=fe_fit1
    fe_template_UV.xarr=sp.xarr
    sp.data=sp2.data - continuous_UV
    
    fe_fit_blue,params1=rescale(sp,"SBB",i,xmin_fe,xbreak,fe_template_UV,mag_order,plot_path=plot_objpath,slim=[0.85,1.0], vlim=[-0,1,0.1], do_fit=fit_fe_UV)
    argbreak=np.argmin(np.abs(sp.xarr.value-xbreak))
    argbalmer=np.argmin(np.abs(sp.xarr.value-xbalmer))
    
    
    fe_fit[:argbreak]=fe_fit_blue[:argbreak]
    fe_fit[argbreak:argbalmer]=fe_fit_red[argbreak:argbalmer]
    fe_fit[argbalmer:]=0.0*fe_fit_red[argbalmer:]
    
    if plot_best_fits:
        copy1=sp.copy() 
        plot_file=plot_objpath + line_to_fit + "_"+ fileroot.split('.')[0] + "_fe_fit.png"

        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        #argxmin=np.argmin(np.abs(sp.xarr.value-xminc))
        #argxmax=np.argmin(np.abs(sp.xarr.value-xmaxc))
        #copy1.crop(argxmin,argxmax)
        pylab.figure()
        #pylab.ylim(ymin=0,ymax=1.1*copy1.data.max())
        #pylab.xlim(xmin=xminc-100.0,xmax=xmaxc+100.0)
 

        #pylab.xlim(xmin=xmin_fe,xmax=xmax_fe)
        #pylab.ylabel(r'$L_{\lambda}\left[10^{'+str(mag_order)+'}$ erg s$^{-1}$  $\AA^{-1}\right]$',fontsize=33)
        pylab.xlabel(r'$\lambda\left[\AA\right]$',fontsize=33)
        pylab.plot(sp.xarr.value,sp2.data-continuous_UV,'k')
        pylab.plot(sp.xarr.value,fe_fit,'r')
        pylab.savefig(plot_file,format='png', dpi=600,bbox_inches='tight')





            
    #-------------Fitting the Blue Bump------------------------
    object_dictionary['MgII_complex']['fwhm']=FWHMV_Mg
    sp.data=sp2.data-continuous_UV-fe_fit
    MgII_fit,object_dictionary['MgII_complex']=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order,plot_path=plot_objpath, linenames=lines_dict['MgII_complex'], do_fit=fit_Mg)
    
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
    sp.xarr.units='angstroms'
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
        
        argxmin=np.argmin(np.abs(sp.xarr.value-xminc))
        argxmax=np.argmin(np.abs(sp.xarr.value-xmaxc))
        copy1.crop(argxmin,argxmax)
    
        pylab.figure()
        pylab.ylim(ymin=0,ymax=1.1*copy1.data.max())
        
        pylab.xlim(xmin=xminc-100.0,xmax=xmaxc+100.0)
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
        
    
    np.savetxt(obj_path +'MgII_fit.txt'         ,np.transpose([sp.xarr.value,MgII_fit]) )
    np.savetxt(obj_path +'feUV_fit.txt'          ,np.transpose([sp.xarr.value,fe_fit]) )
    np.savetxt(obj_path +'continuous_MgII.txt' ,np.transpose([sp.xarr.value,continuous]))
    
    
    
    #np.savetxt(obj_path +'H_fit.txt'     ,np.transpose([sp.xarr.value,H_fit]) )
    
    sample_dictionary[fileroot]['model']    
    

    sample_dictionary[fileroot]['model']['feUV']={}
    

    sample_dictionary[fileroot]['model']['MgII_complex']['datafile']=      obj_path +'MgII_fit.txt'
    sample_dictionary[fileroot]['model']['feUV']['datafile']=              obj_path +'feUV_fit.txt'
    sample_dictionary[fileroot]['model']['feUV']['fwhm']=                  fe_uv_params[0]
    sample_dictionary[fileroot]['model']['feUV']['shift']=                 fe_uv_params[1]#    sample_dictionary[fileroot]['model']['feUV']['scale']=                 fe_uv_params[2]
    sample_dictionary[fileroot]['model']['feUV']['datafile']=              obj_path +'feUV_fit.txt'
    sample_dictionary[fileroot]['model']['MgII_complex']['continuous']={}
    sample_dictionary[fileroot]['model']['MgII_complex']['continuous']['datafile']=        obj_path + 'continuous_MgII.txt'
    sample_dictionary[fileroot]['model']['MgII_complex']['chi2']=        chi2
    sample_dictionary[fileroot]['model']['MgII_complex']['dof']=        dof

    


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
