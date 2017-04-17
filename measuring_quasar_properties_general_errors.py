#to run: python measuring_quasar_properties_general.py line_to_fit path_to_data file_location_with_list_of_spectra_name
#from fitcode.spectral_modeling_lib import *
import copy
import sys
from scipy.integrate import simps as integral
from scipy.stats import kurtosis
from scipy.stats import skew

#matplotlib.use('Agg') 
#matplotlib.use('Agg')

import astropy.units as un
import numpy as np
import pyspeckit
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.interpolate import UnivariateSpline as Interp
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from astropy.cosmology import FlatLambdaCDM
from fitcode.line import *
import json
from matplotlib import rcParams
from scipy.stats import spearmanr as spearman
from scipy.stats import pearsonr as pearson

iterations_per_object=100 # number of interations to calculate errors
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


    
lines_to_fit=[line_to_fit]
savespec=data_dir
dataout_path=data_dir+'fit/'
plots_path= data_dir+'plots/'

filen=dataout_path+'property_dictionary.json'

#---- Typical wavelengths of almost emission-line-free emission where continuum is calculated---#
continuum_bands={}
continuum_bands['SiOIV']=1350.0
continuum_bands['CIV']=1450.0
continuum_bands['MgII']=3000.0
continuum_bands['CIII']=2000.0
continuum_bands['Hbeta']=5100.0
continuum_bands['Halpha']=6200.0
#---- Typical wavelengths of almost emission-line-free emission where continuum is calculated---#


#----------- CONSTANTS TO CALCULATE MBH TAKEN FROM MEJIA-RESTREPO et al 2016 ---------#
#-------logMbh = logK + alpha*log(Llambda/1e44ergs) + slope*log(FWHM/1000kms)---------#
#----------- you can change them for whichever other calibrations------------#
logK={}

logK['CIV']=6.353
logK['MgII']=6.925 
logK['Hbeta']=6.740 
logK['Halpha']=6.891 
logK['SiOIV']=0
logK['CIII']=0

alpha={}

alpha['CIV']=0.599
alpha['MgII']=0.609
alpha['Hbeta']=0.650
alpha['Halpha']=0.634
alpha['SiOIV']=0
alpha['CIII']=0

slope={}

slope['CIV']=2.0
slope['MgII']=2.0
slope['Hbeta']=2.0
slope['Halpha']=2.0
slope['SiOIV']=0
slope['CIII']=0
#-------logMbh = logK + alpha*log(Llambda/1e44ergs) + slope*log(FWHM/1000kms)---------#
#----------- you can change them for whichever other calibrations------------#



#----------- CONSTANTS TO CONVERT  Llambda INTO AN ESTIMATION OF  L5100 TAKEN FROM MEJIA-RESTREPO et al 2016 ---------#
#-------logL5100pred=log(A) + B*log(Llambda). This will assist to calculate Mdot---------#
#---------------- you can change them for whichever other calibrations-------------------#
A={}


A['CIV']=0.56
A['MgII']=0.67 
A['Hbeta']=1.0
A['Halpha']=1.23
A['SiOIV']=0
A['CIII']=0

B={}

B['CIV']=0.88
B['MgII']=0.92 
B['Hbeta']=1.0
B['Halpha']=0.98
B['SiOIV']=0
B['CIII']=0
#-------logL5100pred=log(A) + B*log(Llambda)---------#
#----------- you can change them for whichever other calibrations------------#







object_dictionary={} # Initializing dictionary


"""
To load the dictionary  after you run de code you do the following:

import json
with open(dataout_path+'property_dictionary.json', 'rb') as fp:
    object_dictionary=json.load(fp)
    del(fp)
"""


properties=['FWHM','FWnarrow','FWbroad','SIGMA','EW','DLAMBDA','DLnarrow','DLbroad','DLmax','DL50','DL90','DL95','DLt','Lcont','Lpeak','Lline','Lnarrow','Lbroad','SKEWNESS','KURTOSIS','MASS','MDOT']
values=['mean','median','error_up','error_low','best_fit']

'''
Once it is loaded you can extract any property you want as follows:

object_dictionary['property']['line']['value']
for example:
FWHMCIVmedian=object_dictionary['FWHM']['CIV']['median']
'''


for property in properties:
    object_dictionary[property]={}
    for line_to_fit in lines_to_fit:
        object_dictionary[property][line_to_fit]={}
        for value in values:
            object_dictionary[property][line_to_fit][value]=[]







       
    



incompletes=[]


try:
    sp_to_run=range(len(filenames))
except:
    filenames=np.array([str(filenames)])
    print 'ONLY ONE FILE'
    sp_to_run=[0]



for obj in sp_to_run:


    fileroot= filenames[obj]
        


    obj_path=dataout_path+fileroot + '/'
    try:
        with open(obj_path + line_to_fit +'_dictionary.json', 'rb') as fp:
            obj_dictionary=json.load(fp)
                
    except:
        print obj , " is incomplete"
        incompletes=np.append(incompletes,obj)
        break


        

for obj in sp_to_run:
    if obj in incompletes:
        continue
    
    fileroot= filenames[obj]

    obj_path=dataout_path+fileroot + '/'
    #print fileroot
    #for line_to_fit in ['CIV','MgII','SiOIV','CIII']:

    with open(obj_path + line_to_fit +'_dictionary.json', 'rb') as fp:
        obj_dictionary=json.load(fp)
        """
        try:
        if line_to_fit=='CIV' and obj_dictionary['res_ratio']>0.07:
        break
        except:
        pass
        #print 'object file has been opened ', obj_path + line_to_fit +'_dictionary.json'
        """

       

    spectrum_file=savespec+fileroot#+'.txt'
    sp = pyspeckit.Spectrum(spectrum_file)
        
            
        
       
    continuous_file=savespec+'fit/'+fileroot+'/continuous_'+line_to_fit+'.txt'
    
    xarr,continuous=np.genfromtxt(continuous_file,unpack=1)
        


    amin=np.argmin(np.abs(sp.xarr.value-xarr[0]))
    amax=np.argmin(np.abs(sp.xarr.value-xarr[-1]))
    sp.crop(amin,amax)
    try:
        mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    except:
        print 'problems with the object ', fileroot, ' shoul be related to a problem in the spectrum file itself, not in this code'
        print 'this code will return 0 values for every single quantity related to this object'
        
        for property in properties:
            for line_to_fit in lines_to_fit:
                for value in values:
                    object_dictionary[property][line_to_fit][value]+=[0]
                        

            
                        
        continue
        
            
        
    line_info=global_dict[line_to_fit+'_complex']['lines'][line_to_fit].copy()
    #print line_info
    try:
        model=obj_dictionary['lines']
    except:
        model=obj_dictionary
    continuous=np.interp(sp.xarr,xarr,continuous)
    id=0
    if np.isnan(continuous[id]) or np.isinf(continuous[id]):
        #print 'changing ', id
        continuous[id]=continuous[id+1]
    model1=copy.deepcopy(model)
    #fwhm_array, luminosity_array,EW_array,dlambda_array,std_array,dlmax_array,dl50_array,dl90_array,dl95_array,dlt_array,ske_array,kurto_array,lcont_array,lmax_array,Mass_array,Mdot_array=(np.zeros(iterations_per_object),)*16

    fwhm_array=np.zeros(iterations_per_object)
    luminosity_array=np.zeros(iterations_per_object)
    EW_array=np.zeros(iterations_per_object)
    dlambda_array=np.zeros(iterations_per_object)
    std_array=np.zeros(iterations_per_object)
    dlmax_array=np.zeros(iterations_per_object)
    dl50_array=np.zeros(iterations_per_object)
    dl90_array=np.zeros(iterations_per_object)
    dl95_array=np.zeros(iterations_per_object)
    dlt_array=np.zeros(iterations_per_object)
    ske_array=np.zeros(iterations_per_object)
    kurto_array=np.zeros(iterations_per_object)
    lcont_array=np.zeros(iterations_per_object)
    lmax_array=np.zeros(iterations_per_object)
    Mass_array=np.zeros(iterations_per_object)
    Mdot_array=np.zeros(iterations_per_object)


    
    if True:
        fwhm,fwhm_low,fwhm_up, luminosity,EW,EW_low,EW_up,dlambda,dlambda_low,dlambda_up,lambda0,c_total,conti,varr,xarr,std,ske,kurto,dlmax,dl50,dl90,dl95,dlt=line_parameters(sp,continuous,model1,line_info,linename=line_to_fit)
        fwmin,fwmax,l1,l2,dv1,dv2,fwhmine,fwhmaxe,l1e,l2e,dv1e,dv2e=line_decomposition_measurements(sp,continuous,model,line_info,linename=line_to_fit)
        
        ske=skew(c_total,nan_policy='omit')
        kurto=kurtosis(c_total,nan_policy='omit')
        lmax=c_total.max()*10**(mag_order)*1.0
        
        wavelenght_cont=continuum_bands[line_to_fit]
        arg1450=np.argmin(np.abs(xarr-wavelenght_cont))
        arg14501=np.argmin(np.abs(sp.xarr.value-wavelenght_cont))
        #cont1450=continuous[arg1450]*10**(mag_order)*1450.0*1.0
        cont1450=conti[arg1450]*10**(mag_order)*wavelenght_cont*1.0
        lcont=cont1450
        if arg14501>=10:
            lconte=np.median((sp.error/sp.data)[arg14501-10:arg14501+10])
        else:
            lconte=np.median((sp.error/sp.data)[0:20])

        try:
            lcont_array=np.random.normal(lcont,lconte,iterations_per_object)    
        except:
            lcont_array=lcont*np.ones_like(fwhm_array)
        if line_to_fit=='CIV' or line_to_fit=='MgII':
            print 'MgII or CIV\n '
            fwhm1,fwhm_low1,fwhm_up1, luminosity,EW,EW_low1,EW_up1,dlambda1,dlambda_low1,dlambda_up1,lambda01,c_total1,conti1,varr1,xarr1,std1,ske1,kurto1,dmax1,dl501,dl901,dl951,dlt1=line_parameters(sp,continuous,model,line_info)
            
            lmax=c_total1.max()*10**(mag_order)*1.0
                

                
            
            
        Mass=(logK[line_to_fit]+alpha[line_to_fit]*np.log10(lcont/1e44)+slope[line_to_fit]*np.log10(fwhm/1e3))
        L5100pred=1e44*A[line_to_fit]*(lcont/1e44)**B[line_to_fit]
        Lv5100pred=5100e-8*L5100pred/3e10
        f0=1.2e30# erg/sec/Hz
        bv=2.0
        fth1=f0*0.86*(1+bv*0.86)/(1+bv) #assuming an inclination of 30*
        Mdot=(Lv5100pred/(fth1))**1.5/10**(Mass-8)
        

            
        
        # monte carlo iteration to obtain median and errors
        for itera in range(iterations_per_object):
            for component in line_info['components']:
                if component in model.keys():
                    try:
                        amplitude=np.random.normal(model[component]['modelpars'][0],model[component]['modelerrs'][0])
                    except:
                         amplitude=model[component]['modelpars'][0]
                         
                    try:
                        lambda_center=np.random.normal(model[component]['modelpars'][1],model[component]['modelerrs'][1])
                    except:
                        lambda_center=model[component]['modelpars'][1]
                    try:
                        standard_dev=np.random.normal(model[component]['modelpars'][2],model[component]['modelerrs'][2])
                    except:
                        standard_dev=model[component]['modelpars'][2]
                    model1[component]['modelpars']=[amplitude, lambda_center, standard_dev]
                else:
                    print component, " does not belong to ", line_to_fit, ' complex in the object ', fileroot
                    print "check that you are using the right directory of  the files need for this measurement"
                    print "and or that all the files exists for this particular object"
                    sys.exit()
            fwhm_array[itera], luminosity_array[itera],EW_array[itera],dlambda_array[itera],std_array[itera],dlmax_array[itera],dl50_array[itera],dl90_array[itera],dl95_array[itera],dlt_array[itera]=line_measurements(sp,continuous,model1,line_info,linename=line_to_fit)
            ske_array[itera]=skew(c_total,nan_policy='omit')
            kurto_array[itera]=kurtosis(c_total,nan_policy='omit')
            lmax_array[itera]=c_total.max()*10**(mag_order)*1.0

            
            
            Mass_array[itera]=(logK[line_to_fit]+alpha[line_to_fit]*np.log10(lcont_array[itera]/1e44)+slope[line_to_fit]*np.log10(fwhm_array[itera]/1e3))
            L5100pred=1e44*A[line_to_fit]*(lcont_array[itera]/1e44)**B[line_to_fit]
            Lv5100pred=5100e-8*L5100pred/3e10
            f0=1.2e30# erg/sec/Hz
            bv=2.0
            fth1=f0*0.86*(1+bv*0.86)/(1+bv) #assuming an inclination of 30*
            Mdot_array[itera]=(Lv5100pred/(fth1))**1.5/10**(Mass_array[itera]-8)
            

            


            if line_to_fit=='CIV' or line_to_fit=='MgII':
                fwhm1,fwhm_low1,fwhm_up1, luminosity1,EW1,EW_low1,EW_up1,dlambda1,dlambda_low1,dlambda_up1,lambda01,c_total1,conti1,varr1,xarr1,std1,ske1,kurto1,dmax1,dl501,dl901,dl951,dlt1=line_parameters(sp,continuous,model1,line_info)
                luminosity_array[itera]=luminosity1
                EW_array[itera]=EW1
                lmax_array[itera]=c_total1.max()*10**(mag_order)*1.0

            

    
        bestfit_measurements=np.array([fwhm,fwmin,fwmax,std,EW,dlambda,dv1,dv2,dlmax,dl50,dl90,dl95,dlt,lcont,lmax,luminosity,l1,l2,ske,kurto,Mass,Mdot])#,'Mass','Mdot','LLedd'])
        
        mean_measurements=np.array([np.mean(fwhm_array),fwmin,fwmax,np.mean(std_array),np.mean(EW_array),np.mean(dlambda_array),dv1,dv2,np.mean(dlmax_array),np.mean(dl50_array),np.mean(dl90_array),np.mean(dl95_array),np.mean(dlt_array),np.mean(lcont_array),np.mean(lmax_array),np.mean(luminosity_array),l1,l2,np.mean(ske_array),np.mean(kurto_array),np.mean(Mass_array),np.mean(Mdot_array)])#,'Mass','Mdot','LLedd'])
        
        errorup_measurements=np.array([np.percentile(fwhm_array,84),fwmin+fwhmine,fwmax+fwhmaxe,np.percentile(std_array,84),np.percentile(EW_array,84),np.percentile(dlambda_array,84),dv1e+dv1,dv2e+dv2,np.percentile(dlmax_array,84),np.percentile(dl50_array,84),np.percentile(dl90_array,84),np.percentile(dl95_array,84),np.percentile(dlt_array,84),lcont+lconte,np.percentile(lmax_array,84),np.percentile(luminosity_array,84),l1+l1e,l2+l2e,np.percentile(ske_array,84),np.percentile(kurto_array,84),np.percentile(Mass_array,84),np.percentile(Mdot_array,84)])#,'Mass','Mdot','LLedd'])
        
        errorlow_measurements=np.array([np.percentile(fwhm_array,16),fwmin-fwhmine,fwmax-fwhmaxe,np.percentile(std_array,16),np.percentile(EW_array,16),np.percentile(dlambda_array,16),dv1-dv1e,dv2-dv2e,np.percentile(dlmax_array,16),np.percentile(dl50_array,16),np.percentile(dl90_array,16),np.percentile(dl95_array,16),np.percentile(dlt_array,16),lcont-lconte,np.percentile(lmax_array,16),np.percentile(luminosity_array,16),l1-l1e,l2-l2e,np.percentile(ske_array,16),np.percentile(kurto_array,16),np.percentile(Mass_array,16),np.percentile(Mdot_array,16)])#,'Mass','Mdot','LLedd'])
        
        median_measurements=np.array([np.percentile(fwhm_array,50),fwmin,fwmax,np.percentile(std_array,50),np.percentile(EW_array,50),np.percentile(dlambda_array,50),dv1,dv2,np.percentile(dlmax_array,50),np.percentile(dl50_array,50),np.percentile(dl90_array,50),np.percentile(dl95_array,50),np.percentile(dlt_array,50),lcont,np.percentile(lmax_array,50),np.percentile(luminosity_array,50),l1,l2,np.percentile(ske_array,50),np.percentile(kurto_array,50),np.percentile(Mass_array,50),np.percentile(Mdot_array,50)])#,'Mass','Mdot','LLedd'])#,'Mass','Mdot','LLedd'])
        
        errorup_measurements=errorup_measurements-median_measurements
        errorlow_measurements=median_measurements-errorlow_measurements
        
    if False:
        
        print 'problems with the object ', fileroot, ' shoul be related to a problem in the fit, not in this code'
        print 'this code will return 0 values for every single quantity related to this object'


            
        for property in properties:
            for line_to_fit in lines_to_fit:
                for value in values:
                    object_dictionary[property][line_to_fit][value]+=[0]
                        

        continue
        
        

        
    property_counter=0
    for property in properties:
            
        for line_to_fit in lines_to_fit:
            
            object_dictionary[property][line_to_fit]['mean']+=[mean_measurements[property_counter]]
            object_dictionary[property][line_to_fit]['median']+=[median_measurements[property_counter]]
            object_dictionary[property][line_to_fit]['error_up']+=[errorup_measurements[property_counter]]
            object_dictionary[property][line_to_fit]['error_low']+=[errorlow_measurements[property_counter]]
            object_dictionary[property][line_to_fit]['best_fit']+=[bestfit_measurements[property_counter]]
        property_counter+=1
        

"""    
    if (obj)%1000==0:
        print "loop=  ", obj, ' fileroot ', fileroot
        print "Saving data for safety"
            
            
        if os.path.isfile(filen):
            with open(filen, 'rb') as fp:
                object_dictionary=json.load(fp)
                del(fp)
                
        print "dmax","dl50","dl90","dl95","dlt"
        try:
            print dmax,dl50,dl90,dl95,dlt
            print "FWHM(km/s)=",fwhm        
        except:
            dmax,dl50,dl90,dl95,dlt=0,0,0,0,0
            print dmax,dl50,dl90,dl95,dlt
            print "FWHM(km/s)=",0
        


        print "Saving for Safety. Is the next FWHM array long enough and with proper units??"
        #print FWHMs[line_to_fit]
            
        
        with  open(filen, 'wb') as fp:
            json.dump(object_dictionary, fp)
"""



        


with  open(dataout_path+'property_dictionary_total_'+line_to_fit+'.json', 'wb') as fp:
    json.dump(object_dictionary, fp)

sys.exit()







mhblo=0.65
mhalo=0.45#0.54
mcivlo=0.53
mmgiilo=0.62
bhblo=6+np.log10(5.26)
#bhal=7+np.log10(0.89*2.51)
bhalo=7.30
bcivlo=6.66
bmgiilo=6+np.log10(5.6)


logKc=6.333
alphac=0.599


logKm=6.925
alpham=0.609



logKc=6.66
alphac=0.53

#alpham=0.62
#logKm=6+np.log10(5.6)
logKm=6.925
alpham=0.609



FW=np.array(FWHMs[line_to_fit]) #FWHM(kms)
EW=np.array(EWs[line_to_fit]) # EW equivalent with
Lc=np.array(LCs[line_to_fit]) # luminosity continuum (erg/s)  Lambda*L_lambda
Ll=np.array(Ls[line_to_fit]) # luminosity of the line (erg/s)
Dl=np.array(DLAMBDAs[line_to_fit]) #lineshift centroid (km/s)
Dlm=np.array(DLmax[line_to_fit]) # lineshift of the peak (km/s) -----> DVmax

FWna=np.array(FWn[line_to_fit]) # FWHM of the narrowest component of the BLR
FWbr=np.array(FWb[line_to_fit]) # FWHM of the bradest component of the BLR
Lna=np.array(Ln[line_to_fit]) #  line luminosity of the narrowest component of the BLR
Lbr=np.array(Lb[line_to_fit]) #  line luminosity of the broadest component of the BLR
DLna=np.array(DLn[line_to_fit])  # line shift  of the narrowest component of the BLR
DLbr=np.array(DLb[line_to_fit]) #  line shift of the broadest component of the BLR
DL50max=np.array(DL50[line_to_fit]) # DV50 - DVmax where DV50 is the averaged shift of the central 50 percent of the profile
DL90max=np.array(DL90[line_to_fit]) # DV90 - DVmax where DV90 is the averaged shift of the central 90 percent of the profile
DL95max=np.array(DL95[line_to_fit])# DV95 - DVmax where DV95 is the averaged shift of the central 95 percent of the profile
DLtmax=np.array(DLt[line_to_fit])# DV100 - DVmax where DV100 is the averaged shift of the central 100 percent of the profile


#np.savetxt(dataout_path+'plira_measurements_'+line_to_fit+'.txt',np.transpose([Ll,EW,FW,Dl,Dlm]),header='#Lline(erg/s) EW(AA) FWHM(km/s) DVcentroid(km/s) DVmax(km/s)')
#np.savetxt(dataout_path+'plira_additonal_measurements'+line_to_fit+'.txt',np.transpose([FWna,FWbr,Lna,Lbr,DLna,DLbr,DL50max,DL90max,DL95max,DLtmax]),header='#FWHMnarrow,FWHMbroad,Lnarrow,Lbroad,DLnarrow,DLbroad,DL50percent-DLmax,DL90percent-DLmax,DL95percent-DLmax,DLtotal-DLmax')
if line_to_fit=='CIV':
    FWCIV=np.array(FWHMs[line_to_fit])
    L1450=np.array(CONT1450s[line_to_fit])
    MassCIV=(logKc+alphac*np.log10(L1450/1e44)+2*np.log10(FWCIV/1e3))
    L5100pred=1e44*10**(-0.303)*(L1450/1e44)**0.922
    Lv5100pred=5100e-8*L5100pred/3e10
    f0=1.2e30# erg/sec/Hz
    bv=2.0
    fth1=f0*0.86*(1+bv*0.86)/(1+bv)
    Mdot_realCIV=(Lv5100pred/(fth1))**1.5/10**(MassCIV-8)
    np.savetxt(dataout_path+'MASSES_FWHM_Llambda_CIV.txt',np.transpose([MassCIV,FWCIV,L1450,Mdot_realCIV,L5100pred]),header='#Mass(CIV)[log(Msun)],FWHM(CIV)[km/s],L1450[erg/s]')

if line_to_fit=='MgII':
    FWCIV=np.array(FWHMs[line_to_fit])
    L1450=np.array(CONT3000s[line_to_fit])
    MassCIV=(logKm+alpham*np.log10(L1450/1e44)+2*np.log10(FWCIV/1e3))
    L5100pred=1e44*10**(-0.303)*(L1450/1e44)**0.922
    Lv5100pred=5100e-8*L5100pred/3e10
    f0=1.2e30# erg/sec/Hz
    bv=2.0
    fth1=f0*0.86*(1+bv*0.86)/(1+bv)
    Mdot_realCIV=(Lv5100pred/(fth1))**1.5/10**(MassCIV-8)
    np.savetxt(dataout_path+'MASSES_FWHM_Llambda_MgII.txt',np.transpose([MassCIV,FWCIV,L1450,Mdot_realCIV,L5100pred]),header='#Mass(CIV)[log(Msun)],FWHM(MgII)[km/s],L3000[erg/s]')


