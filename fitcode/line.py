import matplotlib
matplotlib.use('Agg')
import astropy.units as un
from matplotlib import pylab
import numpy as np 
import pyspeckit
import os
import sys
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from scipy.integrate import simps as integral
#import uncertainties as unc  
#import uncertainties.unumpy as unumpy  
from pyspeckit.spectrum.models.template import template_fitter
from scipy import arange, array, exp
from scipy.stats import kurtosis
from scipy.stats import skew


def extrap1d(interpolator):
        xs = interpolator.x
        ys = interpolator.y

        def pointwise(x):
            if x < xs[0]:
                return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            elif x > xs[-1]:
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
                return interpolator(x)

        def ufunclike(xs):
                return array(map(pointwise, array(xs)))
            
        return ufunclike
    


#import quasar_fit_continuum
c=2.99792458e5
def gaussian(x,A0,lambda_0,sigma_l):
    """ definition of the guassian function"""
    if(A0==0 or sigma_l==0):
        y=np.zeros_like(x)
    else:
        y = A0  * np.exp( -( (x - lambda_0 + 0.0) / ( np.sqrt(2.0)*sigma_l ) )**2 )
    return y

def ugaussian(x,A0,lambda_0,sigma_l):
    """ definition of the guassian function"""
    if(A0==0 or sigma_l==0):
        y=unumpy.uarray(( np.zeros_like(x), np.zeros_like(x)))
    else:
        y = A0  * unumpy.exp( -( (x - lambda_0 + 0.0) / ( np.sqrt(2.0)*sigma_l ) )**2 )
    return y

# The following two routines are ony useful for no noisy data
def arglocalmin(data):
    argmin = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1 # local min
    return argmin
def arglocalmax(data):
    argmax = (np.diff(np.sign(np.diff(data))) < 0).nonzero()[0] + 1 # local max
    return argmax
def fwhm_four_comps(sp,continuous,*args):
    #print "before gauss"
    total=0
    s_0=np.empty(3)
    s_up=np.empty(3)
    s_low=np.empty(3)
    luminosity=0
    fwhm=0
    luminosity=0
    for s in args:
        #print len(s)
        #print s
        for k in range(len(s['modelpars'])):
            s_0[k]=s['modelpars'][k]
            
            s_up[k]=s['modelpars'][k] + s['modelerrs'][k]
            
            s_low[k]=s['modelpars'][k] - s['modelerrs'][k]
            luminosity=luminosity+s['flux']
        
        total=total+gaussian(xarr,s_0[0],s_0[1],s_0[2])
        total_up=total+gaussian(sp.xarr,s_up[0],s_up[1],s_up[2])
        total_low=total+gaussian(sp.xarr,s_low[0],s_low[1],s_low[2])
    argmaxi=arglocalmax(total)
    xo=np.linspace(sp.xarr[0],sp.xarr[-1],10000)
    f_obs_func = interp1d(sp.xarr, total, kind='linear')
    if (len(argmaxi)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(sp.xarr, total-total[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xo, f_obs_func(xo))       
        roots_obs=spline_obs.roots()
        #if(len(roots_obs))==0:
        #    fwhm=0
        #else:
        fwhm=roots_obs[-1]-roots_obs[0]
    else:    
        width_obs=0
        reorder=np.argsort(total[argmaxi])
        totalmax=total[argmaxi][reorder][::-1]
        argmaxi=argmaxi[reorder][::-1]
        count=0
        for arg in argmaxi:
            if total[arg]>=0.7*totalmax[0]:
                f_obs_func = interp1d(sp.xarr, total-total[arg]/2.0, kind='linear')
                spline_obs = interpolation(xo, f_obs_func(xo))                            
                roots_obs=spline_obs.roots()
                #if len(roots_obs)==0:
                #    width_obs=width_obs
                #else:
                width_obs=width_obs + roots_obs[-1]-roots_obs[0]
                #print roots_obs[-1],roots_obs[0], totalmax[0]/2
                count=count+1
            
            fwhm=width_obs/count
            
    """
    pylab.figure()
    pylab.plot(sp.xarr,total)
    for m in args:
        ald=gaussian(sp.xarr,m[0],m[1],m[2])
        print "xsize=",np.size(sp.xarr),m
        print "msize",np.size(ald),m
        pylab.plot( sp.xarr, ald  )
            
    pylab.plot( sp.xarr, sp.data-continuous  )
    pylab.plot(roots_obs[0],totalmax[0]/2.0,'ro',markersize=10.0)
    pylab.plot(roots_obs[-1],totalmax[0]/2.0,'ro',markersize=10.0)
    pylab.show()
    """
    
    #luminosity1=integral(spline_obs(xo),x=xo)
    return fwhm, luminosity,total












def line_parameters(sp,continuous,model,line_info,linename='doesnotmatter'):
    """
    Input:
    sp ---> Pyspeckit like spectrum 
    continunuous ---> Fit to the continuum (array must be the same lenght than sp.xarr
    model --> Dictionary of the line complex to work.
    line_info --> It has the info of the lines in the complex. Defined in the configuration file constraints.cfg
    linename --> If the line to measure is a doublet like CIV and MgII please specify the name, otherwise it 'doesnotmatter'
    """
    c=2.99792458e5
    xarr1=np.arange(1100.0,sp.xarr[0].value,(sp.xarr[1].value-sp.xarr[0].value))
    xarr2=np.arange(sp.xarr[-1].value,20000.0,(sp.xarr[1].value-sp.xarr[0].value))
    xarr=np.append(xarr1,sp.xarr.value)
    xarr=np.append(xarr,xarr2)


    #xarrb=(1000.0,10000.0,(sp.xarr[1].value-sp.xarr[0].value)/10.0)
    #cont_obs_func=interp1d(sp.xarr, continuous, kind='linear')
    #cont_obs_func_ext = extrap1d(cont_obs_func)
    #continuous=cont_obs_func_ext(xarr)
    #if continuous[0]==np.nan:
    #continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[1])
    continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[0])
    #print "before gauss"
    total=0.0*xarr
    
    total_up=0
    total_low=0
    s_0=np.empty(3)
    s_up=np.empty(3)
    s_low=np.empty(3)
    luminosity=0
    luminosity_up=0
    luminosity_low=0
    fwhm=0
    fwhm_up=0
    fwhm_low=0
    
    components=line_info['components']
    guesses=line_info['guesses']
    if linename=='CIV' or linename=='MgII':
        components=components[:2]
        guesses=guesses[:2]
        #print 'excluded last two components',components
    #print line_info['components']
    #print model[components[0]]
    #print model[components[1]]
    #exit(0)
    
    lambda0=np.mean(guesses[1::3])
    xarr1=np.arange(1100.0,sp.xarr[0].value,(sp.xarr[1].value-sp.xarr[0].value))
    xarrb=np.arange(lambda0-1000.0,lambda0+1000,0.1)
    totalb=0.0*xarrb
    continuousb=np.interp(xarrb,xarr,continuous,right=continuous[-1],left=continuous[0])
    for component in components:
        #print len(s)
        #print s
        #print component
        if component in model.keys():
            s=model[component]
        else:
            print "pailander"
            return 0,0,0,0,0,0,0,0,0,0,0,0*xarr,0*xarr,0*xarr,0*xarr,0,0,0,0,0,0,0,0

        for k in range(len(s['modelpars'])):
            s_0[k]=s['modelpars'][k]
            
            s_up[k]=s['modelpars'][k] + s['modelerrs'][k]
            if s_up[k]<0:
                s_up[k]=0
		s_low[k]=s['modelpars'][k] - s['modelerrs'][k]
            if s_low[k]<0:
                s_low[k]=0
                s_up[k]=0
            #print s_0    
        luminosity=luminosity+s['flux']
        #luminosity_up=luminosity_up+s_up['flux']
        #luminosity_low=luminosity_low+s_low['flux']
        
        total=total+gaussian(xarr,s_0[0],s_0[1],s_0[2])
	totalb=totalb+gaussian(xarrb,s_0[0],s_0[1],s_0[2])
        #print total
        #exit(0)
        total_up=total_up+gaussian(xarr,s_up[0],s_up[1],s_up[2])
        
        total_low=total_low+gaussian(xarr,s_low[0],s_low[1],s_low[2])
    argmaxi=arglocalmax(totalb)
    argmaxi_up=arglocalmax(total_up)
    argmaxi_low=arglocalmax(total_low)
    xo=np.linspace(xarr[0],xarr[-1],30000)
    f_obs_func = interp1d(xarrb, totalb, kind='linear')
    xb=xarrb
    
    if (len(argmaxi)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarrb, totalb-totalb[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xb, f_obs_func(xb))       
        roots_obs=spline_obs.roots()
        try:
		fwhm=roots_obs[-1]-roots_obs[0]
	except:
		print 'exception 0'
	        return 0,0,0,0,0,0,0,0,0,0,0,0*xarr,0*xarr,0*xarr,0*xarr,0,0,0,0,0,0,0,0
    elif (len(argmaxi)==0 or np.mean(totalb)<=0.0):
        fwhm=0.0
        print "argmax == 0", len(argmaxi)
        #print "mean   max",np.mean(total),"    " ,np.max(total)
    else:    
        width_obs=0
        reorder=np.argsort(totalb[argmaxi])
        totalmax=totalb[argmaxi][reorder][::-1]
        
        argmaxi=argmaxi[reorder][::-1]
        count=0
        #print "argmaxi",argmaxi
        #print np.mean(total)
        
        
        for arg in argmaxi:
            #print "argmax dif= 1", len(argmaxi)
            if totalb[arg]>=0.7*totalmax[0]:
                f_obs_func = interp1d(xarrb, totalb-totalb[arg]/2.0, kind='linear')
                #f_obs_func = interp1d(xarr, total-totalmax[0]/2.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))                  
                
                
                roots_obs=spline_obs.roots()
                roots_obs=np.sort(roots_obs)
                #if len(roots_obs)==0:
                #    width_obs=width_obs
                #else:
                try:
                    width_obs=width_obs + np.abs(roots_obs[-1]-roots_obs[0])
                    count=count+1
                except:
                    width_obs=width_obs
                
        try:
            fwhm=width_obs/(1.0*count)
        except:
            fwhm=width_obs
    if (len(argmaxi_up)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarr, total-total[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xo, f_obs_func(xo))       
        roots_obs=spline_obs.roots()
        #if len(roots_obs)==0:
        #    fwhm_up=0
        #else:
        fwhm_up=roots_obs[-1]-roots_obs[0]
    if (len(argmaxi_up)==0 or np.mean(total_up)<=0.0):
        fwhm_up=0.0
    
    else:    
        width_obs_up=0
        reorder=np.argsort(total_up[argmaxi_up])
        totalmax_up=total_up[argmaxi_up][reorder][::-1]
        
        argmaxi_up=argmaxi_up[reorder][::-1]
        
        count_up=0
    
    
    
              
        for arg_up in argmaxi_up:
            if total_up[arg_up]>=0.7*totalmax_up[0]:
                f_obs_func_up = interp1d(xarr, total_up-total_up[arg_up]/2.0, kind='linear')
                spline_obs_up = interpolation(xo, f_obs_func_up(xo))            
                roots_obs_up=spline_obs_up.roots()
                #if len(roots_obs_up==0):
                #    width_obs_up=width_obs_up
                #else:
                try:
                    width_obs_up=width_obs_up + np.abs(roots_obs[-1]-roots_obs[0])
                    count=count+1
                except:
                    width_obs_up=width_obs_up
                
		
                #width_obs_up=width_obs_up + roots_obs_up[-1]-roots_obs_up[0]
                #count_up=count_up+1
	if count_up==0:
		count_up=1
        fwhm_up=width_obs_up/count_up
    if (len(argmaxi_low)==1):
        #print "argmax = 1"
        try:
            f_obs_func = interp1d(xarr, total-total[argmaxi[0]]/2.0, kind='linear')
	    spline_obs = interpolation(xo, f_obs_func(xo))       
	    roots_obs=spline_obs.roots()
            fwhm_low=roots_obs[-1]-roots_obs[0]
        except:
            fwhm_low=0
    if (len(argmaxi_low)==0 or np.mean(total_low)<=0.0  ):
        fwhm_low=0.0
    
    else:    
        width_obs_low=0
        reorder=np.argsort(total_low[argmaxi_low])
        totalmax_low=total_low[argmaxi_low][reorder][::-1]
        
        argmaxi_low=argmaxi_low[reorder][::-1]
        
        count_low=0

        for arg_low in argmaxi_low:
            if total_low[arg_low]>=0.7*totalmax_low[0]:
                f_obs_func_low = interp1d(xarr, total_low-total_low[arg_low]/2.0, kind='linear')
                spline_obs_low = interpolation(xo, f_obs_func_low(xo))                            
                roots_obs_low=spline_obs_low.roots()
                if len(roots_obs_low)==0:
                    width_obs_low=0
                    count_low=count_low+1
                    continue
                #print "hola",np.mean(total_low)
                #print "roots_obs_low",roots_obs_low

                try:
                    width_obs_low=width_obs_low + np.abs(roots_obs[-1]-roots_obs[0])
                    count=count+1
                except:
                    width_obs_low=width_obs_low
                
                
                #width_obs_low=width_obs_low + roots_obs_low[-1]-roots_obs_low[0]
                #print roots_obs[-1],roots_obs[0], totalmax[0]/2
                #count_low=count_low+1
	if count_low==0:
		count_low=1
        fwhm_low=width_obs_low/count_low
    if fwhm!=0 and (fwhm_low==0):
        fwhm_low=2*fwhm -fwhm_up
        
    EW=integral(totalb/continuousb,xarrb)
    if np.isnan(EW):
            EW=np.sum(totalb*(xarrb[1]-xarrb[0])/continuousb)
    #print total, len(total), np.mean(total), np.mean(total)
    #print continuous, len(continuous), np.mean(continuous)
    #print total/continuous, np.mean(total/continuous)
    #print EW
    EW_up=integral(total_up/continuous,xarr)
    EW_low=integral(total_low/continuous,xarr)
    if EW!=0 and (EW_low==0):
        EW_low=2*EW -EW_up
     
    if np.max(total)<=0 :
        print 'total is a 0 array'
        dlambda=0
        dmax=0
        dl95=0
        dlt=0
        dl50=0
        dl90=0
    else:
        dlambda=(np.sum(totalb*xarrb)/np.sum(totalb) - lambda0)*c/lambda0
        amax=np.argmax(totalb)
        dmax=(xarrb[amax]-lambda0)*c/lambda0
        xb=xarrb    
        f_obs_func = interp1d(xarrb, totalb-totalb[amax]/2.0, kind='linear')
        spline_obs = interpolation(xb, f_obs_func(xb))       
        roots_obs=spline_obs.roots()
        if len(roots_obs)<2:
                dlambda=0
                dmax=0
                dl95=0
                dlt=0
                dl50=0
                dl90=0
        else:
                
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                dl50=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl50=dl50-dmax
                f_obs_func = interp1d(xarrb, totalb-totalb[amax]/10.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))       
                roots_obs=spline_obs.roots()
        
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                
                dl90=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl90=dl90-dmax
                
                f_obs_func = interp1d(xarrb, totalb-totalb[amax]/20.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))       
                roots_obs=spline_obs.roots()
                
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                
                dl95=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl95=dl95-dmax
                dlt=dlambda-dmax

        
        
    if np.max(total_up)<=0 or np.sum(total_up)==0 :
        dlambda_up=0
    else:
        dlambda_up=(np.sum(total_up*xarr)/np.sum(total_up) - lambda0)*c/lambda0
    
    if np.max(total_low)<=0 or np.sum(total_low)==0:
        dlambda_low=0
    else:
        dlambda_low=(np.sum(total_low*xarr)/np.sum(total_low) - lambda0)*c/lambda0
    

    
    if (dlambda!=0 and dlambda_low==0):
        dlambda_low=2*dlambda - dlambda_up
    if (dlambda!=0 and dlambda_up==0):
        dlambda_up=2*dlambda - dlambda_low
    
    #print component,"dlambda=",dlambda_low,dlambda, dlambda_up
    #print "sum",np.sum(total_low) 
    fwhm=fwhm*c/lambda0
    fwhm_up=fwhm_up*c/lambda0
    fwhm_low=fwhm_low*c/lambda0
    xm=np.dot(xarr,total)/np.sum(total)
    x2m=np.dot(xarr*xarr,total)/np.sum(total)
    x3m=np.dot(  (xarr-xm)*(xarr-xm)*(xarr-xm) , total)/(np.sum(total))
    x4m=np.dot(  (xarr-xm)*(xarr-xm)*(xarr-xm)*(xarr-xm) , total)/np.sum(total)
    std1=np.sqrt(np.abs(x2m - xm*xm))
    std=np.sqrt(np.abs(x2m - xm*xm))*c/lambda0
    ske=x3m/np.power(std1,3)
    kurto=x4m/np.power(std1,4) -3 
    varr=(xarr-lambda0)*c/lambda0
    
    
    #luminosity1=integral(spline_obs(xo),x=xo)
    return fwhm,fwhm_low,fwhm_up, luminosity,EW,EW_low,EW_up,dlambda,dlambda_low,dlambda_up,lambda0,total,continuous,varr,xarr,std,ske,kurto,dmax,dl50,dl90,dl95,dlt







def line_measurements(sp,continuous,model,line_info,linename='doesnotmatter'):
    """
    Input:
    sp ---> Pyspeckit like spectrum 
    continunuous ---> Fit to the continuum (array must be the same lenght than sp.xarr
    model --> Dictionary of the line complex to work.
    line_info --> It has the info of the lines in the complex. Defined in the configuration file constraints.cfg
    linename --> If the line to measure is a doublet like CIV and MgII please specify the name, otherwise it 'doesnotmatter'
    """
    c=2.99792458e5
    xarr1=np.arange(1100.0,sp.xarr[0].value,(sp.xarr[1].value-sp.xarr[0].value))
    xarr2=np.arange(sp.xarr[-1].value,20000.0,(sp.xarr[1].value-sp.xarr[0].value))
    xarr=np.append(xarr1,sp.xarr.value)
    xarr=np.append(xarr,xarr2)


    #xarrb=(1000.0,10000.0,(sp.xarr[1].value-sp.xarr[0].value)/10.0)
    #cont_obs_func=interp1d(sp.xarr, continuous, kind='linear')
    #cont_obs_func_ext = extrap1d(cont_obs_func)
    #continuous=cont_obs_func_ext(xarr)
    #if continuous[0]==np.nan:
    #continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[1])
    continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[0])
    #print "before gauss"
    total=0.0*xarr
    
    total_up=0
    total_low=0
    s_0=np.empty(3)
    s_up=np.empty(3)
    s_low=np.empty(3)
    luminosity=0
    luminosity_up=0
    luminosity_low=0
    fwhm=0
    fwhm_up=0
    fwhm_low=0
    
    components=line_info['components']
    guesses=line_info['guesses']
    if linename=='CIV' or linename=='MgII':
        components=components[:2]
        guesses=guesses[:2]
        #print 'excluded last two components',components
    #print line_info['components']
    #print model[components[0]]
    #print model[components[1]]
    #exit(0)
    
    lambda0=np.mean(guesses[1::3])
    xarr1=np.arange(1100.0,sp.xarr[0].value,(sp.xarr[1].value-sp.xarr[0].value))
    xarrb=np.arange(lambda0-1000.0,lambda0+1000,0.1)
    totalb=0.0*xarrb
    continuousb=np.interp(xarrb,xarr,continuous,right=continuous[-1],left=continuous[0])
    for component in components:
        #print len(s)
        #print s
        #print component
        if component in model.keys():
            s=model[component]
        else:
            print "pailander"
	    return 0,0,0,0,0,0,0,0,0,0,0
        for k in range(len(s['modelpars'])):
            s_0[k]=s['modelpars'][k]
            
            s_up[k]=s['modelpars'][k] + s['modelerrs'][k]
            if s_up[k]<0:
                s_up[k]=0
		s_low[k]=s['modelpars'][k] - s['modelerrs'][k]
            if s_low[k]<0:
                s_low[k]=0
                s_up[k]=0
            #print s_0    
        luminosity=luminosity+s['flux']
        #luminosity_up=luminosity_up+s_up['flux']
        #luminosity_low=luminosity_low+s_low['flux']
        
        total=total+gaussian(xarr,s_0[0],s_0[1],s_0[2])
	totalb=totalb+gaussian(xarrb,s_0[0],s_0[1],s_0[2])
        #print total
        #exit(0)
        total_up=total_up+gaussian(xarr,s_up[0],s_up[1],s_up[2])
        
        total_low=total_low+gaussian(xarr,s_low[0],s_low[1],s_low[2])
    argmaxi=arglocalmax(totalb)
    argmaxi_up=arglocalmax(total_up)
    argmaxi_low=arglocalmax(total_low)
    xo=np.linspace(xarr[0],xarr[-1],30000)
    f_obs_func = interp1d(xarrb, totalb, kind='linear')
    xb=xarrb
    
    if (len(argmaxi)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarrb, totalb-totalb[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xb, f_obs_func(xb))       
        roots_obs=spline_obs.roots()
        try:
		fwhm=roots_obs[-1]-roots_obs[0]
	except:
		print 'exception 0'
	        return 0,0,0,0,0,0,0,0,0,0,0
    elif (len(argmaxi)==0 or np.mean(totalb)<=0.0):
        fwhm=0.0
        print "argmax == 0", len(argmaxi)
        #print "mean   max",np.mean(total),"    " ,np.max(total)
    else:    
        width_obs=0
        reorder=np.argsort(totalb[argmaxi])
        totalmax=totalb[argmaxi][reorder][::-1]
        
        argmaxi=argmaxi[reorder][::-1]
        count=0
        #print "argmaxi",argmaxi
        #print np.mean(total)
        
        
        for arg in argmaxi:
            #print "argmax dif= 1", len(argmaxi)
            if totalb[arg]>=0.7*totalmax[0]:
                f_obs_func = interp1d(xarrb, totalb-totalb[arg]/2.0, kind='linear')
                #f_obs_func = interp1d(xarr, total-totalmax[0]/2.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))                  
                
                
                roots_obs=spline_obs.roots()
                roots_obs=np.sort(roots_obs)
                #if len(roots_obs)==0:
                #    width_obs=width_obs
                #else:
                try:
                    width_obs=width_obs + np.abs(roots_obs[-1]-roots_obs[0])
                    count=count+1
                except:
                    width_obs=width_obs
                
        try:
            fwhm=width_obs/(1.0*count)
        except:
            fwhm=width_obs
        
    EW=integral(totalb/continuousb,xarrb)
    if np.isnan(EW):
            EW=np.sum(totalb*(xarrb[1]-xarrb[0])/continuousb)
    #print total, len(total), np.mean(total), np.mean(total)
    #print continuous, len(continuous), np.mean(continuous)
    #print total/continuous, np.mean(total/continuous)
    #print EW 
    if np.max(total)<=0 :
        print 'total is a 0 array'
        dlambda=0
        dmax=0
        dl95=0
        dlt=0
        dl50=0
        dl90=0
    else:
        dlambda=(np.sum(totalb*xarrb)/np.sum(totalb) - lambda0)*c/lambda0
        amax=np.argmax(totalb)
        dmax=(xarrb[amax]-lambda0)*c/lambda0
        xb=xarrb    
        f_obs_func = interp1d(xarrb, totalb-totalb[amax]/2.0, kind='linear')
        spline_obs = interpolation(xb, f_obs_func(xb))       
        roots_obs=spline_obs.roots()
        if len(roots_obs)<2:
                dlambda=0
                dmax=0
                dl95=0
                dlt=0
                dl50=0
                dl90=0
        else:
                
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                dl50=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl50=dl50-dmax
                f_obs_func = interp1d(xarrb, totalb-totalb[amax]/10.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))       
                roots_obs=spline_obs.roots()
        
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                
                dl90=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl90=dl90-dmax
                
                f_obs_func = interp1d(xarrb, totalb-totalb[amax]/20.0, kind='linear')
                spline_obs = interpolation(xb, f_obs_func(xb))       
                roots_obs=spline_obs.roots()
                
                arghm1=np.argmin(np.abs(xarrb-roots_obs[0]))
                arghm2=np.argmin(np.abs(xarrb-roots_obs[-1]))
                
                dl95=(np.sum((totalb*xarrb)[arghm1:arghm2+1])/np.sum(totalb[arghm1:arghm2+1]) - lambda0)*c/lambda0
                dl95=dl95-dmax
                dlt=dlambda-dmax

        
    

    
    #print component,"dlambda=",dlambda_low,dlambda, dlambda_up
    #print "sum",np.sum(total_low) 
    fwhm=fwhm*c/lambda0
    xm=np.dot(xarr,total)/np.sum(total)
    #Higher momentums
    x2m=np.dot(xarr*xarr,total)/np.sum(total)
    #x3m=np.dot(  (xarr-xm)*(xarr-xm)*(xarr-xm) , total)/(np.sum(total))
    #x4m=np.dot(  (xarr-xm)*(xarr-xm)*(xarr-xm)*(xarr-xm) , total)/np.sum(total)
    #std1=np.sqrt(np.abs(x2m - xm*xm))
    std=np.sqrt(np.abs(x2m - xm*xm))*c/lambda0
    #ske=x3m/np.power(std1,3)
    #kurto=x4m/np.power(std1,4) -3 
    varr=(xarr-lambda0)*c/lambda0
    
    
    #luminosity1=integral(spline_obs(xo),x=xo)
    return fwhm, luminosity,EW,dlambda,std,dmax,dl50,dl90,dl95,dlt

















def line_decomposition(sp,continuous,model,line_info,linename='doesnotmatter'):
    """
    Input:
    sp ---> Pyspeckit like spectrum 
    continunuous ---> Fit to the continuum (array must be the same lenght than sp.xarr
    model --> Dictionary of the line complex to work.
    line_info --> It has the info of the lines in the complex. Defined in the configuration file constraints.cfg
    linename --> If the line to measure is a doublet like CIV and MgII please specify the name, otherwise it 'doesnotmatter'
    Output:
    returns the FWHMs, Luminosities and line shifts of the Broad-narrrow (BN)
    and Broad-Broad (BB) components
    fwmin--> FWHM BN
    fwmax--> FWHM BB
    l1 --> Luminosity BN
    l2 --> Luminosity BB
    dv1 --> Offset BN
    dv2 --> Offset BB
    """
    c=2.99792458e5
    xarr1=np.arange(1100.0,sp.xarr[0].value,sp.xarr[1].value-sp.xarr[0].value)
    xarr2=np.arange(sp.xarr[-1].value,20000.0,sp.xarr[1].value-sp.xarr[0].value)
    xarr=np.append(xarr1,sp.xarr.value)
    xarr=np.append(xarr,xarr2)
    
    continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[0])
    #print "before gauss"
    
    s_0=np.empty(3)
    
    luminosity=0
    
    fwhm=0
    
    
    components=line_info['components']
    guesses=line_info['guesses']
    if linename=='CIV' or linename=='MgII':
        components=components[:2]
        guesses=guesses[:2]
        #print 'excluded last two components',components
    
    lambda0=np.mean(guesses[1::3])
    FWHM=[]
    L=[]
    DV=[]
    for component in components:
        #print len(s)
        #print s
        #print component
        if component in model.keys():
            s=model[component]
        else:
            print "pailander"
            return 0,0,0,0,0,0

        
	s_0=np.array(s['modelpars'])
	s_e=np.array(s['modelerrs'])
	luminosity=luminosity+s['flux']

    argfwmin=np.argmin(FWHM)
    argfwmax=np.argmax(FWHM)
    fwhmin=FWHM[argfwmin]
    fwhmax=FWHM[argfwmax]
    L1=L[argfwmin]
    L2=L[argfwmax]
    dv1=DV[argfwmin]
    dv2=DV[argfwmax]
    return fwhmin,fwhmax,L1,L2,dv1,dv2




def line_decomposition_measurements(sp,continuous,model,line_info,linename='doesnotmatter'):
    """
    Input:
    sp ---> Pyspeckit like spectrum 
    continunuous ---> Fit to the continuum (array must be the same lenght than sp.xarr
    model --> Dictionary of the line complex to work.
    line_info --> It has the info of the lines in the complex. Defined in the configuration file constraints.cfg
    linename --> If the line to measure is a doublet like CIV and MgII please specify the name, otherwise it 'doesnotmatter'
    Output:
    returns the FWHMs, Luminosities and line shifts of the Broad-narrrow (BN)
    and Broad-Broad (BB) components
    fwmin--> FWHM BN
    fwmax--> FWHM BB
    l1 --> Luminosity BN
    l2 --> Luminosity BB
    dv1 --> Offset BN
    dv2 --> Offset BB
    """
    c=2.99792458e5
    xarr1=np.arange(1100.0,sp.xarr[0].value,sp.xarr[1].value-sp.xarr[0].value)
    xarr2=np.arange(sp.xarr[-1].value,20000.0,sp.xarr[1].value-sp.xarr[0].value)
    xarr=np.append(xarr1,sp.xarr.value)
    xarr=np.append(xarr,xarr2)
    
    continuous=np.interp(xarr,sp.xarr.value,continuous,right=continuous[-1],left=continuous[0])
    #print "before gauss"
    
    s_0=np.empty(3)
    
    luminosity=0
    
    fwhm=0
    
    
    components=line_info['components']
    guesses=line_info['guesses']
    if linename=='CIV' or linename=='MgII':
        components=components[:2]
        guesses=guesses[:2]
        #print 'excluded last two components',components
    
    lambda0=np.mean(guesses[1::3])
    FWHM=[]
    L=[]
    DV=[]
    FWHMe=[]
    Le=[]
    DVe=[]
    for component in components:
        #print len(s)
        #print s
        #print component
        if component in model.keys():
            s=model[component]
        else:
            print "pailander"
            return 0,0,0,0,0,0

        
	s_0=np.array(s['modelpars'])
	s_e=np.array(s['modelerrs'])
	luminosity=luminosity+s['flux']
	
        FWHM=np.append(FWHM,2.35842*s_0[2]*c/lambda0)
	FWHMe=np.append(FWHMe,2.35842*s_e[2]*c/lambda0)
        DV=np.append(DV,(s_0[1]-lambda0)*c/lambda0)
	DVe=np.append(DVe,(s_e[1])*c/lambda0)
        L=np.append(L,luminosity)
	
	Le=np.append(Le,luminosity*s_e[0]/s_0[0])
	
    argfwmin=np.argmin(FWHM)
    argfwmax=np.argmax(FWHM)
    fwhmin=FWHM[argfwmin]
    fwhmine=FWHMe[argfwmin]
    fwhmax=FWHM[argfwmax]
    fwhmaxe=FWHMe[argfwmax]
    L1=L[argfwmin]
    L1e=Le[argfwmin]
    L2=L[argfwmax]
    L2e=Le[argfwmax]
    dv1=DV[argfwmin]
    dv1e=DVe[argfwmin]
    dv2=DV[argfwmax]
    dv2e=DVe[argfwmax]
    return fwhmin,fwhmax,L1,L2,dv1,dv2,fwhmine,fwhmaxe,L1e,L2e,dv1e,dv2e


















def line_parameters_1(sp,continuous,model,line_info):
    c=2.99792458e5
    xarr1=np.arange(1100.0,sp.xarr[0],sp.xarr[1]-sp.xarr[0])
    xarr=np.append(xarr1,sp.xarr)
    continuous=np.interp(xarr,sp.xarr,continuous)
    #print "before gauss"
    total=0
    total_up=0
    total_low=0
    s_0=np.empty(3)
    s_u=np.empty((3,1000))
    s_up=np.empty(3)
    s_low=np.empty(3)
    luminosity=0
    luminosity_up=0
    luminosity_low=0
    fwhm=0
    fwhm_up=0
    fwhm_low=0
    
    components=line_info['components']
    guesses=line_info['guesses']
    lambda0=np.mean(guesses[1::3])
    for component in components:
        #print len(s)
        #print s
        if component in model.keys():
            s=model[component]
        else:
            return 0,0,0,0,0,0*xarr,0*xarr
        for k in range(len(s['modelpars'])):
            s_0[k]=s['modelpars'][k]
            s_u[k]=np.random.normal(s['modelpars'][k],s['modelerrs'][k],1000)
            s_up[k]=s['modelpars'][k] + s['modelerrs'][k]
            if s_up[k]<0:
                s_up[k]=0
            s_low[k]=s['modelpars'][k] - s['modelerrs'][k]
            if s_low[k]<0:
                s_low[k]=0
            luminosity=luminosity+s['flux']
            #luminosity_up=luminosity_up+s_up['flux']
            #luminosity_low=luminosity_low+s_low['flux']
        
        total=total+gaussian(xarr,s_0[0],s_0[1],s_0[2])
        total_up=total_up+gaussian(xarr,s_up[0],s_up[1],s_up[2])
        
        total_low=total_low+gaussian(xarr,s_low[0],s_low[1],s_low[2])
    argmaxi=arglocalmax(total)
    argmaxi_up=arglocalmax(total_up)
    argmaxi_low=arglocalmax(total_low)
    xo=np.linspace(xarr[0],xarr[-1],10000)
    f_obs_func = interp1d(xarr, total, kind='linear')
    
    if (len(argmaxi)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarr, total-total[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xo, f_obs_func(xo))       
        roots_obs=spline_obs.roots()
        fwhm=roots_obs[-1]-roots_obs[0]
    if (len(argmaxi)==0 or np.mean(total)<=0.0):
        fwhm=0.0
    else:    
        width_obs=0
        reorder=np.argsort(total[argmaxi])
        totalmax=total[argmaxi][reorder][::-1]
        
        argmaxi=argmaxi[reorder][::-1]
        count=0
        #print "argmaxi",argmaxi
        #print np.mean(total)
        
        for arg in argmaxi:
            if total[arg]>=0.7*totalmax[0]:
                f_obs_func = interp1d(xarr, total-total[arg]/2.0, kind='linear')
                spline_obs = interpolation(xo, f_obs_func(xo))                  
                
                
                roots_obs=spline_obs.roots()
                
                
                width_obs=width_obs + roots_obs[-1]-roots_obs[0]
                count=count+1
                
        fwhm=width_obs/(1.0*count)
                
    if (len(argmaxi_up)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarr, total-total[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xo, f_obs_func(xo))       
        roots_obs=spline_obs.roots()
        fwhm_up=roots_obs[-1]-roots_obs[0]
    if (len(argmaxi_up)==0 or np.mean(total_up)<=0.0):
        fwhm_up=0.0
    
    else:    
        width_obs_up=0
        reorder=np.argsort(total_up[argmaxi_up])
        totalmax_up=total_up[argmaxi_up][reorder][::-1]
        
        argmaxi_up=argmaxi_up[reorder][::-1]
        
        count_up=0
    
    
    
              
        for arg_up in argmaxi_up:
            if total_up[arg_up]>=0.7*totalmax_up[0]:
                f_obs_func_up = interp1d(xarr, total_up-total_up[arg_up]/2.0, kind='linear')
                spline_obs_up = interpolation(xo, f_obs_func_up(xo))            
                roots_obs_up=spline_obs_up.roots()
                width_obs_up=width_obs_up + roots_obs_up[-1]-roots_obs_up[0]
                count_up=count_up+1
        fwhm_up=width_obs_up/count_up
    if (len(argmaxi_low)==1):
        #print "argmax = 1"
        f_obs_func = interp1d(xarr, total-total[argmaxi[0]]/2.0, kind='linear')
        spline_obs = interpolation(xo, f_obs_func(xo))       
        roots_obs=spline_obs.roots()
        fwhm_low=roots_obs[-1]-roots_obs[0]
    if (len(argmaxi_low)==0 or np.mean(total_low)<=0.0  ):
        fwhm_low=0.0
    
    else:    
        width_obs_low=0
        reorder=np.argsort(total_low[argmaxi_low])
        totalmax_low=total_low[argmaxi_low][reorder][::-1]
        
        argmaxi_low=argmaxi_low[reorder][::-1]
        
        count_low=0

        for arg_low in argmaxi_low:
            if total_low[arg_low]>=0.7*totalmax_low[0]:
                f_obs_func_low = interp1d(xarr, total_low-total_low[arg_low]/2.0, kind='linear')
                spline_obs_low = interpolation(xo, f_obs_func_low(xo))                            
                roots_obs_low=spline_obs_low.roots()
                if len(roots_obs_low)==0:
                    width_obs_low=0
                    count_low=count_low+1
                    continue
                #print "hola",np.mean(total_low)
                #print "roots_obs_low",roots_obs_low
                
                width_obs_low=width_obs_low + roots_obs_low[-1]-roots_obs_low[0]
                #print roots_obs[-1],roots_obs[0], totalmax[0]/2
                count_low=count_low+1

        fwhm_low=width_obs_low/count_low
    if fwhm!=0 and (fwhm_low==0):
        fwhm_low=2*fwhm -fwhm_up
        
    EW=integral(total/continuous,xarr)
    EW_up=integral(total_up/continuous,xarr)
    EW_low=integral(total_low/continuous,xarr)
    if EW!=0 and (EW_low==0):
        EW_low=2*EW -EW_up
    
    if np.max(total)<=0 :
        dlambda=0
    else:
        dlambda=(np.sum(total*xarr)/np.sum(total) - lambda0)*c/lambda0
    
    
    if np.max(total_up)<=0 or np.sum(total_up)==0 :
        dlambda_up=0
    else:
        dlambda_up=(np.sum(total_up*xarr)/np.sum(total_up) - lambda0)*c/lambda0
    
    if np.max(total_low)<=0 or np.sum(total_low)==0:
        dlambda_low=0
    else:
        dlambda_low=(np.sum(total_low*xarr)/np.sum(total_low) - lambda0)*c/lambda0
    

    
    if (dlambda!=0 and dlambda_low==0):
        dlambda_low=2*dlambda - dlambda_up
    if (dlambda!=0 and dlambda_up==0):
        dlambda_up=2*dlambda - dlambda_low
    
    print component,"dlambda=",dlambda_low,dlambda, dlambda_up
    print "sum",np.sum(total_low) 
    fwhm=fwhm*c/lambda0
    fwhm_up=fwhm_up*c/lambda0
    fwhm_low=fwhm_low*c/lambda0
    varr=(xarr-lambda0)*c/lambda0
    """
    pylab.figure()
    pylab.plot(xarr,total)
    for m in args:
        ald=gaussian(xarr,m[0],m[1],m[2])
        print "xsize=",np.size(xarr),m
        print "msize",np.size(ald),m
        pylab.plot( xarr, ald  )
            
    pylab.plot( xarr, sp.data-continuous  )
    pylab.plot(roots_obs[0],totalmax[0]/2.0,'ro',markersize=10.0)
    pylab.plot(roots_obs[-1],totalmax[0]/2.0,'ro',markersize=10.0)
    pylab.show()
    """
    
    #luminosity1=integral(spline_obs(xo),x=xo)
    return fwhm,fwhm_low,fwhm_up, luminosity,EW,EW_low,EW_up,dlambda,dlambda_low,dlambda_up,lambda0,total,varr



















def fwhm_two_comps(xarr,CIV11par,CIV12par):

    s1=gaussian(xarr,CIV11par[0],CIV11par[1],CIV11par[2])
    s2=gaussian(xarr,CIV12par[0],CIV12par[1],CIV12par[2])
    
    total=s1+s2
    argmax=arglocalmax(total)
    

    xo=np.linspace(xarr[0],xarr[-1],10000)
    f_obs_func = interp1d(xarr, total, kind='cubic')

    if (len(argmax)==1):
        spline_obs = interpolation(xo, f_obs_func(xo)-total[argmax[0]]/2.0, s=0)
        roots=spline_obs.roots()
        fwhm=roots_obs[-1]-roots_obs[0]
    else:
        fwhm=0
        for arg in argmax:
            spline_obs = interpolation(xo, f_obs_func(xo)-total[arg]/2.0, s=0)
            roots=spline_obs.roots()
            width_obs=width_obs + roots_obs[-1]-roots_obs[0]
        fwhm=width_obs/len(argmax)
    return fwhm





def kms_to_wl(kms,line_center):
    c=3e5
    wl=kms*line_center/c
    return wl

def wl_to_kms(wl,line_center):
    c=3e5
    kms=wl*c/line_center
    return kms

def galex_mAB_Fnuv(mAB):
    Fnuv=2.06e-16*np.power(10, (20.08-mAB)/2.5 )
    return Fnuv

def galex_mAB_Ffuv(mAB):
    Ffuv=1.40e-15*np.power(10, (18.82-mAB)/2.5 )
    return Ffuv

    


def continuous_substraction( spec_number, sp, magorder,cont_limits,region):
    
    model=np.loadtxt("./model/cont_model_" + str(spec_number) +".txt")
    
    Lmodel=fdata=np.interp(sp.xarr.value,model[:,0],model[:,1])
    
    Fmodel=np.multiply( Lmodel ,10**(-1*magorder) )
    model[:,1]=np.multiply( model[:,1] ,10**(-1*magorder) )
    
    #printA Lmodel,Fmodel, sp.data, magorder
    
    exclude_file = "./excludes/exclude_cont.txt"
    
    
    exclude_cont=np.loadtxt(exclude_file,skiprows=2)
    
    cont=exclude_cont[spec_number,:][1:]
    cont=cont[cont<8999]
    print "cont=  ", cont
    arg_min=np.argmin( np.abs( cont - cont_limits[0] ) )
    if cont[arg_min] - cont_limits[0] < 0 : arg_min= arg_min + 1
    arg_max=np.argmin( np.abs( cont - cont_limits[1] ) )
    if cont[arg_max] - cont_limits[0] > 0 : arg_max= arg_max - 1
    if len(cont[arg_min:arg_max+1]) % 2 == 1: arg_max=arg_max + 1
    print len(cont[arg_min:arg_max+1]), len(cont[arg_min:arg_max+1]) % 2
    cont=cont[arg_min:arg_max+1]
    
    if region=="OP":
        print "-------------------REGION OP--------------------"
        exclude_file = "./excludes/exclude_cont_OP.txt"
        cont=np.loadtxt(exclude_file,skiprows=2)
        cont=cont[spec_number,:]
        #cont=[4467.0, 4474.0]
    
    if region=="OPnear":
        print "-------------------REGION OP--------------------"
        exclude_file = "./excludes/exclude_cont_OPnear.txt"
        cont=np.loadtxt(exclude_file,skiprows=2)
        cont=cont[spec_number,:]
        #cont=[4467.0, 4474.0]
    if region=="OPfar":
        print "-------------------REGION OPfar--------------------"
        exclude_file = "./excludes/exclude_cont_OPfar.txt"
        cont=np.loadtxt(exclude_file,skiprows=2)
        cont=cont[spec_number,:]
        #cont=[4467.0, 4474.0]
        
    print "cont=  ", cont
    
    

    Nregions=np.int( len(cont)/2 )
    
    print "Nreg= ", Nregions

    arg_list=np.ones(len(cont),dtype="int")
    i=0

    for wl in cont:
        arg= np.argmin(np.abs(sp.xarr.value - wl) ) 
        arg_list[i]=arg
        i=i+1
    ratio=np.zeros(Nregions)
    
    for i in range(Nregions):
        ratio[i]=np.mean(sp.data[ arg_list[2*i]:arg_list[2*i + 1] ]/np.mean(Fmodel[ arg_list[2*i]:arg_list[2*i + 1] ]) )
        
        mean_ratio=np.mean(ratio)
    
    
    Fmodel*=mean_ratio
    print "ratio","mean_ratio"
    print ratio,mean_ratio
    if mean_ratio<0: mean_ratio*=-1
    model[:,1]=model[:,1]*mean_ratio

    
    return Fmodel, cont[0], model


#sp (the spectrum) must be given with the global continuous substracted
def balmer_normalization(sp,balmer_cont_template,balmer_lines_template):
    
    template_file=balmer_cont_template
    balmer_lines_file=balmer_lines_template
    balmer_template = np.loadtxt(template_file)
    balmer_lines=np.loadtxt(balmer_lines_file)
    balmer_data=np.interp(sp.xarr.value,balmer_template[:,0], balmer_template[:,1])
    balmer_lines_data=np.interp(sp.xarr.value,balmer_lines[:,0], balmer_lines[:,1])
    balmer_template=pyspeckit.Spectrum(data=balmer_data,xarr=sp.xarr)
    balmer_template.xarr.unit='angstroms'
    balmer_template.xarr.units='angstroms'
    balmer_template.xarr.xtype='wavelength'
    balmer_lines=pyspeckit.Spectrum(data=balmer_lines_data,xarr=sp.xarr)
    balmer_lines.xarr.unit='angstroms'
    balmer_lines.xarr.units='angstroms'
    balmer_lines.xarr.xtype='wavelength'
    
    lambda_balmer=3640
    index=np.argmin( np.abs( lambda_balmer - sp.xarr) )
    index_hi=np.argmin( np.abs( 4000.0 - sp.xarr) )
    flux=np.mean(sp.data[ index - 30 : index + 30 ])
    balmer_lines.data*=flux/balmer_template.data[-1]
    balmer_template.data*=flux/balmer_template.data[-1]
            
    balmer_tot=np.zeros(len(balmer_template.data))
    balmer_tot[:index]=balmer_template.data[:index]
    balmer_tot[index:index_hi]=balmer_lines.data[index:index_hi]
    
    
    return balmer_template, balmer_tot, index
    
def total_continuous_fit(spec_number,sp,continuous,balmer,magorder,region,plot_path="./plots/"):
    

    
    #-----------------------------defining balmer template fitter ------------------------------------#
    

    
    def total_cont(xarr, scale_cont,scale_balmer,cont=continuous.copy(), balm=balmer.copy()):
            
            
        total_continuous =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                          cont.xarr.as_unit('angstroms'),
                                                          scale_cont*cont.data + scale_balmer*balm.data)
        
                    
        return total_continuous
    
        
    modelclass = pyspeckit.spectrum.models.model.SpectralModel(total_cont,2, 
                                                               parnames=['scale_cont','scale_balmer'],
                                                               parlimited=[(True,True),(True,True)],
                                                               parlimits=[(0.6,1.2),(0.85,1.15)],
                                                               shortvarnames=('Ac','Ab') 
                                                               )
    modelclass.__name__ ="contmodel"
        
    xmin=2150
    xmax=3670
    #-----------restarting plotter------------#
    #sp.plotter.refresh()    
    #sp.plotter.reset_limits()
    #sp.plotter(figure=1,xmin=xmin,xmax=xmax)
    #-----------restarting plotter------------#
    
    #exclude=[1000,2190,2200,3641,3651,10000]
    #exclude=[1000,2147,2152,3641,3651,10000]
    exclude_file = "./excludes/exclude_cont_UV.txt"
    exclude=np.loadtxt(exclude_file,skiprows=2)
    exclude=exclude[spec_number,:]
    
    
    
        
    #sp.Registry.add_fitter('total',total,2,multisingle='milti')
    #sp.specfit(fittype='al',exclude=exclude,guesses=[0.99,0.99],xmin=xmin,xmax=xmax)
    

        
    sp.Registry.add_fitter('contmodel',modelclass,2,multisingle='multi')
    sp.specfit(fittype='contmodel',exclude=exclude,guesses=[1.0,1.0],xmin=xmin,xmax=xmax)
    

    #sp.specfit.plot_components(add_baseline=False)
    #sp.plotter.refresh()
    #sp.specfit.plot_fit(annotate=False)
    
    # -----------------  fitting ---------------------
        
     # -----------------  ploting  fit---------------------

    
    fitmin=sp.specfit.xmin
    fitmax=sp.specfit.xmax
    sp.specfit.fullsizemodel() 
    print sp.specfit.parinfo.values[1], sp.specfit.parinfo.values[0], sp.specfit.xmin, sp.specfit.xmax
    
    #--------------rescaling balmer and continuous----------
    continuous.data*=sp.specfit.parinfo.values[0]
    balmer.data*=sp.specfit.parinfo.values[1]
    scale_balmer=sp.specfit.parinfo.values[1]


    #plot_file=plot_path + "new_cont" + line_name + "_" + str(spec_number) + ".png"
    #plot_file=plot_path + "new_cont.png"
    #pylab.figure()
    #pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
    #pylab.xlim(xmin=xmin-300,xmax=xmax)
    #pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
    #pylab.xlabel(r'$\AA$')
    #pylab.plot(sp.xarr,sp.data,'k')
    #pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
    #pylab.savefig(plot_file)
    #pylab.close()
    # -----------------  ploting fit---------------------


    
    model=sp.specfit.fullmodel

    return continuous, balmer, scale_balmer















###----------------------Quasar Continuum + Host Galaxy -----------------######


def quasar_host_continuous(spec_number,sp,continuous,galactic_fraction,magorder,region,plot_path="./plots/"):
    

    
    #-----------------------------defining balmer template fitter ------------------------------------#
    sp1=sp.copy()
    scale_galaxy=np.arange(galactic_fraction-0.1,galactic_fraction +0.1,0.01)
    chi2=np.empty_like(scale_galaxy)
    scale_cont=np.empty_like(scale_galaxy)
    i=0
    for s_gal in scale_galaxy:
        arglambda=np.argmin(np.abs(continuous.xarr-6100.0))
        sp.data=sp1.data - s_gal*continuous.data[arglambda]*np.ones_like(sp.data)
        
        print sp.data
        
        def total_cont(xarr, scale_cont,cont=continuous.copy()):
            
            arglambda=np.argmin(np.abs(cont.xarr-6100.0))
            total_continuous =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                              cont.xarr.as_unit('angstroms'),
                                                              scale_cont*cont.data)
            
        
        

        
            return total_continuous
    
    
        modelclass = pyspeckit.spectrum.models.model.SpectralModel(total_cont,1, 
                                                               parnames=['scale_cont'],
                                                               parlimited=[(True,True)],
                                                               parlimits=[(0.5,1.2)],
                                                               shortvarnames=('Ac') 
                                                               )
        modelclass.__name__ ="agnhost"
        
        xmin=4000
        xmax=7400
        #-----------restarting plotter------------#
        sp.plotter.refresh()    
        sp.plotter.reset_limits()
        sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        #-----------restarting plotter------------#
    
 
        exclude_file = "./excludes/exclude_host.txt"
        exclude=np.loadtxt(exclude_file,skiprows=1)
        exclude=exclude[spec_number,:]
    
    
        sp.Registry.add_fitter('agnhost',modelclass,1,multisingle='multi')
        sp.specfit(fittype='agnhost',exclude=exclude,guesses=[1.0],xmin=xmin,xmax=xmax)
        #chi2[i]=sp.specfit.chi2
        chi2[i]=sp.specfit.optimal_chi2()#sp.specfit.chi2
        scale_cont[i]=sp.specfit.parinfo.values[0]
        i=i+1

    argmin=np.argmin(chi2)
    chi2min=chi2[argmin]
    scont_min=scale_cont[argmin]
    sgal_min=scale_galaxy[argmin]
    galactic_cont=sgal_min*continuous.data[arglambda]*np.ones_like(sp.data)
    sp.data=sp1.data
    #sp.data=sp1.data - sgal_min*continuous.data[arglambda]*np.ones_like(sp.data) 
    #sp.specfit.plot_components(add_baseline=False)
    #sp.plotter.refresh()
    #sp.specfit.plot_fit(annotate=False)
    print chi2
    # -----------------  fitting ---------------------
        
    # -----------------  ploting  fit---------------------
    
    
    fitmin=sp.specfit.xmin
    fitmax=sp.specfit.xmax
    sp.specfit.fullsizemodel() 
    print scont_min, sgal_min, sp.specfit.xmin, sp.specfit.xmax
    
    #--------------rescaling balmer and continuous----------
    
    
    #sp.specfit.fullmodel=continuous.data
    
    arglambda=np.argmin(np.abs(sp.xarr-6100.0))
    #pylab.yscale("log")
    #plot_file=plot_path + "new_cont" + line_name + "_" + str(spec_number) + ".png"
    plot_file=plot_path + "cont_host.png"
    pylab.figure()
    pylab.ylim(ymin=-1.2*sp.data[fitmin:fitmax].min(),ymax=1.3*sp.data[fitmin:fitmax].max())
    pylab.xlim(xmin=4000,xmax=7400)
    pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr,continuous.data*scont_min,'b')
    pylab.plot(sp.xarr,continuous.data[arglambda]*sgal_min*np.ones_like(sp.xarr),'gray')
    pylab.plot(sp.xarr,continuous.data*scont_min + continuous.data[arglambda]*sgal_min*np.ones_like(continuous.data),'r')
    pylab.savefig(plot_file)
    #pylab.close()
    # -----------------  ploting fit---------------------


    continuous.data=continuous.data*scont_min + continuous.data[arglambda]*sgal_min*np.ones_like(continuous.data)
    model=continuous.data*scont_min + continuous.data[arglambda]*sgal_min*np.ones_like(continuous.data)
    out_pars=np.array([scont_min,sgal_min])
    return continuous,galactic_cont,out_pars







###----------------------Quasar Continuum + Host Galaxy -----------------######



def smooth(data,smooth,smoothtype='gaussian',downsample=False,downsample_factor=None,
        convmode='same'):
    """
    Smooth and downsample the data array

    Parameters
    ----------
    smooth  :  float 
        Number of pixels to smooth by
    smoothtype : [ 'gaussian','hanning', or 'boxcar' ]
        type of smoothing kernel to use
    downsample :  bool 
        Downsample the data?
    downsample_factor  :  int 
        Downsample by the smoothing factor, or something else?
    convmode : [ 'full','valid','same' ]
        see :mod:`numpy.convolve`.  'same' returns an array of the same length as
        'data' (assuming data is larger than the kernel)
    """

    
    roundsmooth = smooth # can only downsample by integers

    if downsample_factor is None and downsample:
        downsample_factor = roundsmooth
    elif downsample_factor is None:
        downsample_factor = 1
    downsample_factor = 1
    if smooth > len(data) or downsample_factor > len(data):
        raise ValueError("Error: trying to smooth by more than the spectral length.")

    
    
    xkern  = np.linspace(-5*smooth,5*smooth,smooth*11)
    kernel = np.exp(-xkern**2/(2*(smooth/np.sqrt(8*np.log(2)))**2))
    kernel /= kernel.sum()
    if len(kernel) > len(data):
        lengthdiff = len(kernel)-len(data)
        if lengthdiff % 2 == 0: # make kernel same size as data
            kernel = kernel[lengthdiff/2:-lengthdiff/2]
        else: # make kernel 1 pixel smaller than data but still symmetric
            kernel = kernel[lengthdiff/2+1:-lengthdiff/2-1]
    
    # deal with NANs or masked values
    if hasattr(data,'mask'):
        if type(data.mask) is np.ndarray:
            OK = True - data.mask
            if OK.sum() > 0:
                data = pyspeckit.interpolation._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])
            else:
                data = OK
    if np.any(True - np.isfinite(data)):
        OK = np.isfinite(data)
        if OK.sum() > 0:
            data = pyspeckit.interpolation._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])
        else:
            data = OK

    if np.any(True - np.isfinite(data)):
        raise ValueError("NANs in data array after they have been forcibly removed.")

    smdata = np.convolve(data,kernel,convmode)#[::downsample_factor]

    return smdata







    
def blue_bump_fitter(sp,line_name,spec_number,xmin,xmax,fe_template, magorder,plot_path="./plots/", do_fit=True):
    if do_fit:
        
        number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        del(number,mag,redshift)
        
    #-----------------------------defining fe template fitter ------------------------------------#
        
        for i in range(3):
        
            def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

                new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                           fe_template.xarr.as_unit('angstroms')+shift_fe, 
                                                                           fe_template.data )
           
                new_template=  new_fe_template #+ new_balmer_template 
                return new_fe_template


        


            modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 

                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
            modelclass.__name__ = "template"



        #-----------------------------defining fe template fitter ------------------------------------#

        
        #-----------restarting plotter------------#
            sp.plotter.refresh()    
            sp.plotter.reset_limits()
            sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        #-----------restarting plotter------------#
        
        
        # -----------------  fitting ---------------------
        
            exclude_file = "./excludes/exclude_" + line_name + ".txt"
            if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
                exclude_data=np.loadtxt(exclude_file,skiprows=2)
                exclude=exclude_data[spec_number,:]
            else:
                exclude=[]
            
            
        
            sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
            sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        
            sp.specfit.fullsizemodel()  
        
            fe_template.data*=sp.specfit.parinfo.values[0]
            fe_template.xarr+=sp.specfit.parinfo.values[1]
            new_template=fe_template.copy()
        
            def suavizado(xarr, smooth_factor_fe,templat=new_template.copy()):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
            
            
                data=smooth(templat.data,smooth_factor_fe , downsample=False)   
                templat.data=data
                n_template =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                            templat.xarr.as_unit('angstroms'), 
                                                            templat.data)
                fe_template.data=data
                    
                return n_template
    
        
            soft = pyspeckit.spectrum.models.model.SpectralModel(suavizado,
                                                                 1, 
                                                                 parnames=['smooth_fe'],
                                                                 parlimited=[(True,False)],
                                                                 parlimits=[(0,0)],
                                                                 shortvarnames=(r'sigma_x') )
            soft.__name__ = "suavizado"
        
            sp.Registry.add_fitter('suavizado',soft,1,multisingle='multi')
            sp.specfit(fittype='suavizado',exclude=exclude,guesses=[1.0],xmin=xmin,xmax=xmax)
        
        # -----------------  fitting ---------------------
        
    # -----------------  ploting  fit---------------------

    
        fitmin=sp.specfit.xmin
        fitmax=sp.specfit.xmax
        sp.specfit.fullsizemodel() 
        
        plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=xmin-300,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr,sp.data,'k')
        pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
        pylab.savefig(plot_file)
        #pylab.close()
        # -----------------  ploting fit---------------------
    
    
       # ----------plotting residuals-------------
        
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
        pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
            
        pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
        plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.savefig(plot_file)
        #pylab.close()
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model

def line_fitter(sp,line_name, spec_number,guesses, limits, limited, tied, xmin, xmax,plot_path="./plots/", offset=-0.4, magorder=51,linenames=[],excluding=0 ,do_fit=True):
    if do_fit:
        
        exclude=[]
        # -----------------  fitting ---------------------
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if excluding:
            if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
		    exclude_data=np.loadtxt(exclude_file,skiprows=2)
		    exclude=exclude_data[spec_number,:]
	    else:
		    exclude=[]

        
        sp.specfit(exclude=exclude,guesses = guesses, limits=limits, limited=limited, tied=tied, annotate = False,fittype='gaussian',xmin=xmin,xmax=xmax)
        #sp.specfit(exclude=exclude,guesses = guesses, limits=limits, limited=limited, tied=tied, annotate = False,fittype='lorentzian',xmin=xmin,xmax=xmax,multifit=False)
        # -----------------  fitting ---------------------
        if False:
             
            #-----------restarting plotter------------#
            sp.plotter.refresh()    
            sp.plotter.reset_limits()
            sp.plotter(figure=1,xmin=xmin,xmax=xmax)
            #-----------restarting plotter------------#

            # ----------------- fit plotting  ---------------------
            sp.specfit.plot_components(add_baseline=False)
            pylab.close()
            sp.plotter.refresh()
            sp.specfit.plot_fit(annotate=False)
            # -----------------  fit plotting  ---------------------

            # -----------------  line measuring  ---------------------
            
            # ----------------- TO DEVELOP!!!!  ---------------------
            # ----------------- SiOIV line measuring  ---------------------

            #----------plotting residuals-------------
            plot_file=plot_path + line_name + "_"  + ".png"
            sp.plotter.figure.savefig(plot_file)
            pylab.close()
            fitmin=np.argmin(np.abs(sp.xarr-xmin))
            fitmax=np.argmin(np.abs(sp.xarr-xmax))
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals.max())
            pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
            pylab.ylabel(r'$10^{'+str(magorder)+'}$ erg s$^{-1}$  $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr[sp.specfit.xmin:sp.specfit.xmax],sp.specfit.residuals,'k')
            plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number] + ".png"
    
            pylab.savefig(plot_file)
            pylab.close()
        # ----------plotting residuals-------------#
    
    
       
    

        # Print out spectral line information
        #filename=plot_path + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".txt"
        #f=open(filename,"w")

        line_dictionary={}
        sp.measure(z = 0.00, fluxnorm = 10**(magorder))
        a=0
        
        for line in sp.measurements.lines.keys():
            
            line_dictionary[linenames[a]]=sp.measurements.lines[line]
            a=a+1
            #f.write("#Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)")
            #print_line = '{0}  {1}  {2}  {3} {4} {5}\n'.format(line, sp.measurements.lines[line]['flux'],sp.measurements.lines[line]['amp'],sp.measurements.lines[line]['fwhm'],sp.measurements.lines[line]['lum'], sp.measurements.lines[line]['modelpars'][1] )
            #f.write(print_line)
        #f.close()
        #---------------measuring----------------------#
       
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
        line_dictionary={}
    return model, line_dictionary
    
    
def fe_fitter(sp,line_name,spec_number,xmin,xmax,fe_template, magorder,plot_path="./plots/", do_fit=True,lambda0=2800.0):
    if do_fit:
        
        #number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        #del(number,mag,redshift)
        
        #-----------------------------defining fe template fitter ------------------------------------#
        
        
        
        def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                
            
            new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.value, 
                                                                           fe_template.xarr.value+shift_fe, 
                                                                           fe_template.data )
                
            new_template=  new_fe_template #+ new_balmer_template 
            return new_fe_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                       
                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
        modelclass.__name__ = "template"

            

        #-----------------------------defining fe template fitter ------------------------------------#

        
        
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
        
        sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
        sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
	
        fe_template.data*=sp.specfit.parinfo.values[0]
        
        wltemp=fe_template.xarr.value+sp.specfit.parinfo.values[1]
        fe_template.xarr=wltemp

	#print fe_template.xarr.unit

        new_template=fe_template.copy()            
        deltav=50.0/2.355
        sigmatemplate=900.0/2.355
            
            
        templates=np.empty( (len(fe_template.data),200) )
        chi2=np.empty(200)
        amplitude=np.empty(200)
        shift=np.empty(200)
        fwhm=np.empty(200)
        
        templates[:,0]=fe_template.data
        
        chi2[0]=sp.specfit.chi2
        amplitude[0]=sp.specfit.parinfo.values[0]
        shift[0]=sp.specfit.parinfo.values[1]
        fwhm[0]=sigmatemplate*2.355
        pos0=np.argmin(np.abs(sp.xarr.value-lambda0))    
        for i in range(1,200):
            
            dl=sp.xarr[pos0+1].value-sp.xarr[pos0].value
            delta_sigma=np.sqrt(2*deltav*sigmatemplate+deltav**2)
            sigma=kms_to_wl(delta_sigma,lambda0)/(2.355*dl)
            
            templates[:,i]=gauss_conv(templates[:,i-1], sigma, order=0, output=None, mode='reflect', cval=0.0)
            fe_template.data=templates[:,i]
            sigmatemplate=sigmatemplate+deltav
            
            def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

                new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
                                                                           fe_template.xarr+shift_fe, 
                                                                           fe_template.data )
                    
                return new_fe_template


        


            modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
            modelclass.__name__ = "template1"


        
            
            
            
            sp.Registry.add_fitter('template1',modelclass,3,multisingle='multi')
            sp.specfit(fittype='template1',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax,debug=False,verbose=False)
            
            sp.specfit.fullsizemodel()  
            
            
            
            chi2[i]=sp.specfit.chi2
            amplitude[i]=sp.specfit.parinfo.values[0]
            shift[i]=sp.specfit.parinfo.values[1]
            fwhm[i]=sigmatemplate*2.355
                
        arg_chi=np.argmin(chi2)
        #arg_chi=0
        chi2_min=chi2[arg_chi]
        amplitude_min=amplitude[arg_chi]
        shift_min=shift[arg_chi]
        fwhm_min=fwhm[arg_chi]
        new_template.data=templates[:,arg_chi]
        new_template.data*=amplitude_min
        new_template.xarr+=shift_min
        
	out_params=np.empty(3)    
        out_params[0]=fwhm_min
        out_params[1]=shift_min
        out_params[2]=amplitude_min
        
        
	print fe_template.unit
	#new_template.xarr.set_unit(un.Angstrom)
	#print sp.xarr.unit
        #new_template = pyspeckit.interpolation.interp(new_template,sp)
	new_template1=sp.copy()
	new_template1.data=np.interp(sp.xarr,new_template.xarr,new_template.data,left=0,right=0)
        new_template=new_template1.copy()
	sp.specfit.residuals=sp.data - new_template.data
        sp.specfit.fullmodel=new_template.data
        
        
        # -----------------  fitting ---------------------
        
    # -----------------  ploting  fit---------------------
        
        sp.specfit.fullsizemodel() 
        fitmin=np.argmin(np.abs(sp.xarr.value-xmin))
        fitmax=np.argmin(np.abs(sp.xarr.value-xmax))
        
        if False:
            #-----------restarting plotter------------#
            sp.plotter.refresh()    
            sp.plotter.reset_limits()
            sp.plotter(figure=1,xmin=xmin,xmax=xmax)
            #-----------restarting plotter------------#
        
        
            plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
            pylab.xlim(xmin=xmin-300,xmax=xmax)
            pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr,sp.data,'k',label="FWHM="+str(fwhm_min))
            pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
            pylab.savefig(plot_file)
            pylab.legend()
            #pylab.close()
            # -----------------  ploting fit---------------------
        
        
            # --- -------plotting residuals-------------
        
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
            pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
            pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
            plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number] + ".png"
            pylab.savefig(plot_file)
            #pylab.close()
            # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model,out_params




def fe_scale(sp,fe_uv_params,line_name,spec_number,xmin,xmax,fe_template, continuous, galaxy_template,fitter_name,magorder,lambda0,objectname='test',plot_path="./plots/"  ,do_fit=True,  plot_figure=0,central_wl=5000.0):
    if do_fit:
        sp_copy=sp.copy()
        sp.plotter.autorefresh=False
        #number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        #del(number,mag,redshift)
        
        new_template=fe_template.copy()            
        new_template1 = pyspeckit.interpolation.interp(new_template,sp)
        deltav=100.0/2.355
        sigmatemplate=900.0/2.355
            
        x2800=np.argmin(np.abs(sp.xarr.value-central_wl))
        dl=sp.xarr.value[x2800+1]-sp.xarr.value[x2800]
        sigma_model=fe_uv_params[0]/2.355
        delta_sigma=np.sqrt(sigma_model*sigma_model-sigmatemplate*sigmatemplate)
        sigma=kms_to_wl(delta_sigma,lambda0)/dl 
        
        fe_template.data=gauss_conv(fe_template.data, sigma, order=0, output=None, mode='reflect', cval=0.0)
        arg_fe=np.argmin(np.abs(sp.xarr.value - 5100.0))
        if spec_number==27:arg_fe=np.argmin(np.abs(sp.xarr.value - 4550.0))
        fe_template.data*=np.abs(np.mean(sp.data[arg_fe-10:arg_fe+10])/np.mean(new_template1.data[arg_fe-10:arg_fe+10]))
        sp_backup=sp.copy()
        sp1=sp.copy()
        sp1.data=sp.data+continuous.data
        
        fe_template.xarr=fe_template.xarr+fe_uv_params[1]*fe_template.xarr.unit
        
        
        
        #-----------------------------defining fe template fitter ------------------------------------#
        #-----------restarting plotter------------#
        
        #sp.plotter.refresh()    
        #sp.plotter.reset_limits()
        #sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        
        #-----------restarting plotter------------#
       
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
            
        else:
            exclude=[]
            
        s_cont=np.arange(1.0,1.01,0.02)
        shift_total=np.arange(0.0,1.0,2.5)
        chi2=np.empty(len(shift_total)*len(s_cont))
        amplitude=np.empty(len(shift_total)*len(s_cont))
        shift=np.empty(len(shift_total)*len(s_cont))
        scale_c=np.empty(len(shift_total)*len(s_cont))
                
        i=0
        
        for shift_fe in shift_total:
            for scale_cont in s_cont:
                fitter_name='scale'+str(np.random.rand(1)[0])
                s_gal=0.0
                fe_template1=fe_template.copy()
                #fe_template1.xarr=fe_template.xarr + shift_fe
                #sp.data=sp1.data-continuous.data*scale_cont
                
                def model(xarr, scale_fe,dl_fe):
                
			try:
				new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
											   fe_template.xarr.value+shift_fe, 
											   fe_template.data )
			except:    
				new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
											   fe_template.xarr+shift_fe, 
											   fe_template.data )
			
                    
			new_template=  new_fe_template 
			return new_fe_template
                

                


                modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                           2, 
                                                                           parnames=['scale_fe','shift_fe'],
                                                                           parlimited=[(True,True),(True,True)],
                                                                           parlimits=[(0.003,10.0),(-16,16.0)],
                                                                           shortvarnames=(r'Af',r'dv'),
                                                                           )
                modelclass.__name__ =fitter_name
            
                #-----------------------------defining fe template fitter ------------------------------------#

                
            
            
        
                sp.Registry.add_fitter(fitter_name,modelclass,2,override=True,multisingle='multi')
                sp.specfit(fittype=fitter_name,exclude=exclude,guesses=[1.0,0.0],xmin=xmin,xmax=xmax,debug=False,verbose=False,interactive=False)
                
                """
                template_fitter = pyspeckit.models.template_fitter(fe_template1,xshift_units='angstroms')
                sp.Registry.add_fitter('template1',template_fitter,2,override=True,multisingle='multi')
                sp.specfit(fittype='template1',guesses=[1.0,0.0],limits=[(0.1,10.0),(-20.0,20.0)],limited=[(True,True),(True,True)],exclude=exclude)

                """


                """
                
                # Create the fitter from the template spectrum and "Register" it
                t_fitter = template_fitter(fe_template1,xshift_units='angstroms')
                sp.Registry.add_fitter('template',t_fitter,2,override=True,multisingle='multi')

                # The fitted parameters are amplitude & xshift
                # perform the fit:
                sp.specfit(fittype='template',guesses=[1.0,0.0])

                # print the results
                print "scale", "shift"
                print sp.specfit.parinfo
                """

                sp.specfit.fullsizemodel()  
                chi2[i]=sp.specfit.chi2
                amplitude[i]=sp.specfit.parinfo.values[0]
                shift[i]=sp.specfit.parinfo.values[0]#shift_fe
                scale_c[i]=scale_cont
                
        
                i=i+1
            
            
            
        
        
        # -----------------  fitting ---------------------
            
        # -----------------  ploting  fit---------------------
        
                
        arg_chi=np.argmin(chi2)
        
        chi2_min=chi2[arg_chi]
        amplitude_min=amplitude[arg_chi]
        shift_min=shift[arg_chi]
        scale_c_min=scale_c[arg_chi]
        #scale_g_min=scale_g[arg_chi]
        
        continuous.data=scale_c_min*continuous.data
        
        sp.data=sp1.data - continuous.data #- galaxy_template.data
        
        fe_template.data=amplitude_min*fe_template.data
        fe_template.xarr=fe_template.xarr+shift_min*fe_template.xarr.unit
        
        out_params=np.empty(3)    
        out_params[0]=amplitude_min
        out_params[1]=shift_min
        out_params[2]=scale_c_min
        
        
        sorted=np.argsort(chi2)

        new_template = pyspeckit.interpolation.interp(fe_template,sp)
        
        sp.specfit.residuals=sp.data - new_template.data
        sp.specfit.fullmodel=new_template.data
                
        
        sp.specfit.fullsizemodel() 
        fitmin=np.argmin(np.abs(sp.xarr.value-xmin))
        fitmax=np.argmin(np.abs(sp.xarr.value-xmax))

        if plot_figure:
            plot_file=plot_path + "fe_fit_" + line_name + "_"  + objectname+".png"
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
            pylab.xlim(xmin=xmin-300,xmax=xmax)
            pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr,sp.data,'k')
            pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
            pylab.savefig(plot_file)
            
            
        
        
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(np.min(sp.data[fitmin:fitmax])),ymax=1.1*np.max(sp.specfit.residuals[fitmin:fitmax]))
            pylab.xlim(xmin=sp.xarr.value[fitmin],xmax=sp.xarr.value[fitmax])
            pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            
            pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
            plot_file=plot_path + line_name +"_res_" +"_"+ objectname + '.png' #pylab.close()
            pylab.savefig(plot_file)
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model,continuous,fe_template,out_params,chi2_min

def fe_galaxy(sp,fe_uv_params,line_name,spec_number,xmin,xmax,fe_template, continuous, galaxy_template,magorder,plot_path="./plots/", do_fit=True):
    if do_fit:
        
        number,name,mag,group,redshift,galactc_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        del(number,mag,redshift)
        
        new_template=fe_template.copy()            
        new_template1 = pyspeckit.interpolation.interp(new_template,sp)
        #deltav=100.0/2.355
        #sigmatemplate=900.0/2.355
            
        
        
        #sigma_model=fe_uv_params[0]/2.355
        #delta_sigma=np.sqrt(sigma_model*sigma_model-sigmatemplate*sigmatemplate)
        #sigma=kms_to_wl(delta_sigma,4680) 
        
        #fe_template.data=gauss_conv(fe_template.data, sigma, order=0, output=None, mode='reflect', cval=0.0)
        #arg_fe=np.argmin(np.abs(sp.xarr - 5200.0))
        #if spec_number==27:arg_fe=np.argmin(np.abs(sp.xarr - 4550.0))
        #fe_template.data*=np.mean(sp.data[arg_fe-10:arg_fe+10])/np.mean(new_template1.data[arg_fe-10:arg_fe+10])
        sp_backup=sp.copy()
        sp1=sp.copy()
        sp1.data=sp1.data+continuous.data
        arg6700=np.argmin(np.abs(continuous.xarr - 6700.0))
        arg6800=np.argmin(np.abs(continuous.xarr - 6800.0))
        galaxy_template.data*=np.mean(continuous.data[arg6700:arg6800])/np.mean(galaxy_template.data[arg6700:arg6800])
        galaxy_template.data=np.mean(continuous.data[arg6700:arg6800])/np.mean(galaxy_template.data[arg6700:arg6800])*np.ones_like(galaxy_template.data)
        fe_template.xarr+=fe_uv_params[1]
        
        
        #-----------------------------defining fe template fitter ------------------------------------#
        #-----------restarting plotter------------#
        
        sp.plotter.refresh()    
        sp.plotter.reset_limits()
        sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        
        #-----------restarting plotter------------#
       
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
        s_cont=np.arange(0.70,1.0,0.02)
        s_gal=np.arange(0.0,3.0,0.1)
        
        
        
        chi2=np.empty(len(s_cont)*len(s_gal))
        amplitude=np.empty(len(s_cont)*len(s_gal))
        shift=np.empty(len(s_cont)*len(s_gal))
        scale_c=np.empty(len(s_cont)*len(s_gal))
        scale_g=np.empty(len(s_cont)*len(s_gal))
        
        i=0
        sp.data=sp1.data #- scale_cont*continuous.data - scale_galaxy*galaxy_template.data
        def model_gal(xarr, scale_fe,scale_cont, scale_galaxy):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
            
            
            new_gal_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('Angstroms'), 
                                                                       fe_template.xarr.as_unit('Angstroms'), 
                                                                       fe_template.data ) + scale_cont*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                                                                                       continuous.xarr.as_unit('angstroms'), 
                                                                                                                                       continuous.data ) + scale_galaxy*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                                                                                                                                                         galaxy_template.xarr.as_unit('angstroms'), 
                                                                                                                                                                                                         galaxy_template.data )
            
                    
            return new_gal_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model_gal,
                                                                   3, 
                                                                   parnames=['scale_fe','scale_cont','scale_galaxy'],#,'scale_balmer'], 
                                                                   parlimited=[(True,True),(True,True),(True,True)],#,(False,False)],  
                                                                   parlimits=[(0.4,2.5), (0.4,1.0),(0.1,3.0)],#,(0.9,1.1)], 
                                                                   shortvarnames=(r'Af',r'Ac',r'Ag'),#,r'Ab')
                                                                    )
        modelclass.__name__ = "gal_template"

                #-----------------------------defining fe template fitter ------------------------------------#
        
        
            
            
        
        sp.Registry.add_fitter('gal_template',modelclass,3,multisingle='multi')
        sp.specfit(fittype='gal_template',exclude=exclude,guesses=[1.0,0.7,0.3],xmin=xmin,xmax=xmax,debug=False,verbose=False)
        
        
        sp.specfit.fullsizemodel()  
        chi2_min=sp.specfit.optimal_chi2()#sp.specfit.chi2
        amplitude_min=sp.specfit.parinfo.values[0]
        scale_c_min=sp.specfit.parinfo.values[1]
        scale_g_min=sp.specfit.parinfo.values[2]
        
                
                
            
            
            
        
        
        # -----------------  fitting ---------------------
            
        # -----------------  ploting  fit---------------------
        
                
        
        
        continuous.data=scale_c_min*continuous.data
        galaxy_template.data=scale_g_min*galaxy_template.data
        sp.data=sp1.data - continuous.data - galaxy_template.data
        
        fe_template.data*=amplitude_min
        
        
        out_params=np.empty(3)    
        out_params[0]=amplitude_min
        out_params[1]=scale_c_min
        out_params[2]=scale_g_min
        print "s_fe   ,   scont  ,  sgal =  ", amplitude_min , "  ,   "   ,  scale_c_min, "  ,  ", scale_g_min
        sorted=np.argsort(chi2)
        
        
        
        new_template = pyspeckit.interpolation.interp(fe_template,sp)

        sp.specfit.residuals=sp.data - new_template.data
        sp.specfit.fullmodel=new_template.data
        
        
        
        
        
        
        sp.specfit.fullsizemodel() 
        fitmin=np.argmin(np.abs(sp.xarr-xmin))
        fitmax=np.argmin(np.abs(sp.xarr-xmax))
        print fitmin, fitmax, len(sp.data)
        
        plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=xmin-300,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr,sp.data,'k')
        model=sp.specfit.fullmodel - galaxy_template.data - continuous.data
        pylab.plot(sp.xarr,model,'r')
        pylab.savefig(plot_file)
        
            
        #pylab.close()
        # -----------------  ploting fit---------------------
        
    
        # ----------plotting residuals-------------
        
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(np.min(sp.data[fitmin:fitmax])),ymax=1.1*np.max(sp.specfit.residuals[fitmin:fitmax]))
        pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
        pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
            
        pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
        plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number]  #pylab.close()
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        
    else:
        model=0.0*sp.data
    return model,continuous, fe_template,galaxy_template,out_params


def fe_cont_galaxy(sp,fe_uv_params,line_name,spec_number,xmin,xmax,fe_template, continuous, galaxy_template,magorder,plot_path="./plots/", do_fit=True):
    if do_fit:
        
        number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        del(number,mag,redshift)
        
        new_template=fe_template.copy()            
        new_template1 = pyspeckit.interpolation.interp(new_template,sp)
        #deltav=100.0/2.355
        #sigmatemplate=900.0/2.355
            
        
        
        #sigma_model=fe_uv_params[0]/2.355
        #delta_sigma=np.sqrt(sigma_model*sigma_model-sigmatemplate*sigmatemplate)
        #sigma=kms_to_wl(delta_sigma,4680) 
        
        #fe_template.data=gauss_conv(fe_template.data, sigma, order=0, output=None, mode='reflect', cval=0.0)
        #arg_fe=np.argmin(np.abs(sp.xarr - 5200.0))
        #if spec_number==27:arg_fe=np.argmin(np.abs(sp.xarr - 4550.0))
        #fe_template.data*=np.mean(sp.data[arg_fe-10:arg_fe+10])/np.mean(new_template1.data[arg_fe-10:arg_fe+10])
        sp_backup=sp.copy()
        sp1=sp.copy()
        sp1.data=sp1.data+continuous.data
        arg6700=np.argmin(np.abs(continuous.xarr - 6700.0))
        arg6800=np.argmin(np.abs(continuous.xarr - 6800.0))
        print galaxy_template.data
        print np.mean(continuous.data[arg6700:arg6800])
        print np.mean(galaxy_template.data[arg6700:arg6800])
        
        galaxy_template.data*=np.mean(continuous.data[arg6700:arg6800])/np.mean(galaxy_template.data[arg6700:arg6800])
        galaxy_template.data=np.mean(continuous.data[arg6700:arg6800])/np.mean(galaxy_template.data[arg6700:arg6800])*np.ones_like(galaxy_template.data)
        print galaxy_template.data
        print np.mean(continuous.data[arg6700:arg6800])
        print np.mean(galaxy_template.data[arg6700:arg6800])
        
        #sys.exit(1)
        #fe_template.xarr+=fe_uv_params[1]
        
        
        #-----------------------------defining fe template fitter ------------------------------------#
        #-----------restarting plotter------------#
            
        sp.plotter.refresh()    
        sp.plotter.reset_limits()
        sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        
        #-----------restarting plotter------------#
       
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
        s_cont=np.arange(0.96,1.0,0.02)
        s_gal=np.arange(0.18,0.22,0.02)
        #s_gal=np.array([0.0])
        #s_cont=[1.0]
        
        chi2=np.empty(len(s_cont)*len(s_gal))
        amplitude=np.empty(len(s_cont)*len(s_gal))
        shift=np.empty(len(s_cont)*len(s_gal))
        scale_c=np.empty(len(s_cont)*len(s_gal))
        scale_g=np.empty(len(s_cont)*len(s_gal))
        
        i=0
        
        for scale_cont in s_cont:
            for scale_galaxy in s_gal:
                #print sp.data[:30]
                #print continuous.data[:30]
                #print galaxy_template.data[:30]
                sp.data=sp1.data - scale_cont*continuous.data - scale_galaxy*galaxy_template.data
                def model1(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

                    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                               fe_template.xarr.as_unit('angstroms'), 
                                                                               fe_template.data )
                    
                    new_template=  new_fe_template #+ new_balmer_template 
                    return new_fe_template


        


                modelclass = pyspeckit.spectrum.models.model.SpectralModel(model1,
                                                                           3, 
                                                                           parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                           parlimited=[(True,True),(True,True),(True,False)],#,(False,False)],  
                                                                           parlimits=[(0.3,1.5), (-2.0,2.0),(0,0)],#,(0.9,1.1)], 
                                                                           shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                           centroid_par='shift_fe' )
                modelclass.__name__ = "template1"
            
                #-----------------------------defining fe template fitter ------------------------------------#

        
            
            
        
                sp.Registry.add_fitter('template1',modelclass,3,multisingle='multi')
                sp.specfit(fittype='template1',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax,debug=False,verbose=False)

                
                sp.specfit.fullsizemodel()  
                chi2[i]=sp.specfit.optimal_chi2()#sp.specfit.chi2
                amplitude[i]=sp.specfit.parinfo.values[0]
                shift[i]=0.0#sp.specfit.parinfo.values[1]
                scale_c[i]=scale_cont
                scale_g[i]=scale_galaxy
        
                i=i+1
            
            
            
        
        
        # -----------------  fitting ---------------------
            
        # -----------------  ploting  fit---------------------
        
                
        arg_chi=np.argmin(chi2)
        sp.specfit.fullsizemodel() 
        chi2_min=chi2[arg_chi]
        amplitude_min=amplitude[arg_chi]
        shift_min=shift[arg_chi]
        scale_c_min=scale_c[arg_chi]
        scale_g_min=scale_g[arg_chi]
        
        continuous.data=scale_c_min*continuous.data
        galaxy_template.data=scale_g_min*galaxy_template.data
        sp.data=sp1.data - continuous.data - galaxy_template.data
        
        fe_template.data*=amplitude_min
        fe_template.xarr+=shift_min
        
        out_params=np.empty(4)    
        out_params[0]=amplitude_min
        out_params[1]=shift_min
        out_params[2]=scale_c_min
        out_params[3]=scale_g_min
        print "sfe,    scont  ,  sgal =  ",amplitude_min,"  ,  "  ,scale_c_min, "  ,  ", scale_g_min
        sorted=np.argsort(chi2)
        print scale_c[sorted][:10]
        print scale_g[sorted][:10]
        print np.sort(chi2)[:10]
        
        new_template = pyspeckit.interpolation.interp(fe_template,sp)

        sp.specfit.residuals=sp.data - new_template.data
        
        
        
        
        
        
        
        
        sp.specfit.fullmodel=new_template.data
        fitmin=np.argmin(np.abs(sp.xarr-xmin))
        fitmax=np.argmin(np.abs(sp.xarr-xmax))
        print fitmin, fitmax, len(sp.data)
        
        plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=xmin-300,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr,sp.data,'k')
        pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
        pylab.savefig(plot_file)
        
            
        #pylab.close()
        # -----------------  ploting fit---------------------
        
    
        # ----------plotting residuals-------------
        
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(np.min(sp.data[fitmin:fitmax])),ymax=1.1*np.max(sp.specfit.residuals[fitmin:fitmax]))
        pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
        pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
            
        pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
        plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number]  #pylab.close()
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model,continuous,fe_template,galaxy_template,out_params













    
def fe_refitter(sp,line_name,spec_number,xmin,xmax,fe_template,mg_broad, magorder,plot_path="./plots/", do_fit=True):
    if do_fit:
        
        number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        del(number,mag,redshift)
        
    #-----------------------------defining fe template fitter ------------------------------------#
        
        
        
        def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

            new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                           fe_template.xarr.as_unit('angstroms')+shift_fe, 
                                                                           fe_template.data )
                
            new_template=  new_fe_template #+ new_balmer_template 
            return new_fe_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                       
                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
        modelclass.__name__ = "template"

            

        #-----------------------------defining fe template fitter ------------------------------------#

        
        #-----------restarting plotter------------#
        sp.plotter.refresh()    
        sp.plotter.reset_limits()
        sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        #-----------restarting plotter------------#
        
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
        
        sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
        sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
        fe_template.data*=sp.specfit.parinfo.values[0]
        fe_template.xarr+=sp.specfit.parinfo.values[1]
        
        new_template=fe_template.copy()            
        deltav=100.0/2.355
        sigmatemplate=900.0/2.355
        deltavo=(mg_broad - 500)/2.355  - sigmatemplate   
        
        templates=np.empty( (len(fe_template.data),100) )
        chi2=np.empty(10)
        amplitude=np.empty(10)
        shift=np.empty(10)
        fwhm=np.empty(10)
        templates[:,0]=fe_template.data
        #chi2[0]=sp.specfit.optimal_chi2()
        chi2[0]=sp.specfit.chi2
        amplitude[0]=sp.specfit.parinfo.values[0]
        shift[0]=sp.specfit.parinfo.values[1]
        fwhm[0]=sigmatemplate*2.355
        
        delta_sigma=np.sqrt(2*deltavo*sigmatemplate+deltavo**2)
        sigma=kms_to_wl(delta_sigma,2800) 
        templates[:,0]=gauss_conv(templates[:,0], sigma, order=0, output=None, mode='reflect', cval=0.0)
        fe_template.data=templates[:,0]
        sigmatemplate=(mg_broad-500)/2.355
        

    
        for i in range(1,10):
                
            delta_sigma=np.sqrt(2*deltav*sigmatemplate+deltav**2)
            sigma=kms_to_wl(delta_sigma,2800) 
            templates[:,i]=gauss_conv(templates[:,i-1], sigma, order=0, output=None, mode='reflect', cval=0.0)
            fe_template.data=templates[:,i]
            sigmatemplate=sigmatemplate+deltav
            
            def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

                new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                           fe_template.xarr.as_unit('angstroms')+shift_fe, 
                                                                           fe_template.data )
                    
                return new_fe_template


        


            modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
            modelclass.__name__ = "template"


        
            
            
        
            sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
            sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax,debug=False,verbose=False)
            
            sp.specfit.fullsizemodel()  
            
            #chi2[i]=sp.specfit.chi2
            
            chi2[i]=sp.specfit.chi2
            amplitude[i]=sp.specfit.parinfo.values[0]
            shift[i]=sp.specfit.parinfo.values[1]
            fwhm[i]=sigmatemplate*2.355
                
        arg_chi=np.argmin(chi2)
        arg_chi=28
        chi2_min=chi2[arg_chi]
        amplitude_min=amplitude[arg_chi]
        shift_min=shift[arg_chi]
        fwhm_min=fwhm[arg_chi]
        new_template.data=templates[:,arg_chi]
        new_template.data*=amplitude_min
        new_template.xarr+=shift_min
        outf_params=np.empty(3)    
        out_params[0]=fwhm_min
        out_params[1]=shift_min
        out_params[2]=amplitude_min
        
        print np.argsort(chi2)[:10]

        new_template = pyspeckit.interpolation.interp(new_template,sp)
                    
        sp.specfit.residuals=sp.data - new_template.data
        sp.specfit.fullmodel=new_template.data
        
        print "arg_chi,chi,amp,shift,fwhm= ",arg_chi,chi2_min,amplitude_min,shift_min,fwhm_min
        # -----------------  fitting ---------------------
        
        # -----------------  ploting  fit---------------------

        sp.specfit.fullsizemodel() 
        fitmin=np.argmin(np.abs(sp.xarr-xmin))
        fitmax=np.argmin(np.abs(sp.xarr-xmax))
        print fitmin, fitmax, len(sp.data)
        
        plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=xmin-300,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr,sp.data,'k',label="FWHM="+str(fwhm_min))
        pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
        pylab.savefig(plot_file)
        pylab.legend()
        #pylab.close()
        # -----------------  ploting fit---------------------
    
    
       # ----------plotting residuals-------------
        
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
        pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')

            
        pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
        plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.savefig(plot_file)
        #pylab.close()
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model,out_params











    
def rescale(sp,line_name,spec_number,xmin,xmax,fe_template, magorder,plot_path="./plots/", slim=[0.1,10.0],vlim=[-10.0,10.0],do_fit=True):
    if do_fit:
        
        #number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        #del(number,mag,redshift)
        
    #-----------------------------defining fe template fitter ------------------------------------#
        
        
        
        def model(xarr, scale_fe,shift_fe):
                
            try:
		    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr.value+shift_fe, 
									       fe_template.data )
            except:    
		        new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr+shift_fe, 
									       fe_template.data )
            new_template=  new_fe_template 
            return new_fe_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                   2, 
                                                                   parnames=['scale_fe','shift_fe'],#,'scale_balmer'], 
                                                                   parlimited=[(True,True),(True,True)],
                                                                   parlimits=[(slim[0],slim[1]), (vlim[0],vlim[1])],
                                                                   
                                                                   shortvarnames=(r'Af',r'\Delta xf'),#,r'Ab'),
                                                                   centroid_par='shift_fe' )
        fitter_name="rescale"+str(np.random.rand(1)[0])
        modelclass.__name__ = fitter_name

            

        #-----------------------------defining fe template fitter ------------------------------------#

        
        
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
            
        sp.Registry.add_fitter(fitter_name,modelclass,2,override=True,multisingle='multi')
        sp.specfit(fittype=fitter_name,exclude=exclude,guesses=[1.0,0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
        fe_template.data*=sp.specfit.parinfo.values[0]
        
	try:
		wltemp=fe_template.xarr.value+sp.specfit.parinfo.values[1]
		fe_template.xarr=wltemp
	except:
		fe_template.xarr+=sp.specfit.parinfo.values[1]
        
        
        
        fitmin=np.argmin(np.abs(sp.xarr.value-xmin))
        fitmax=np.argmin(np.abs(sp.xarr.value-xmax))
        print fitmin, fitmax, len(sp.data)
        
        if False:
            #-----------restarting plotter------------#
            sp.plotter.refresh()    
            sp.plotter.reset_limits()
            sp.plotter(figure=1,xmin=xmin,xmax=xmax)
            #-----------restarting plotter------------#
        
            plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
            pylab.xlim(xmin=xmin-300,xmax=xmax)
            pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr,sp.data,'k')
            pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
            pylab.savefig(plot_file)
            
            #pylab.close()
            # -----------------  ploting fit---------------------
        
        sp.specfit.fullsizemodel() 
        sp.data=sp.data-sp.specfit.fullmodel
    return sp.specfit.fullmodel, sp.specfit.parinfo.values





    
def scale_and_shift(sp,line_name,spec_number,xmin,xmax,fe_template, magorder,plot_path="./plots/", slim=[0.1,10.0],vlim=[-10.0,10.0],do_fit=True):
    if do_fit:
        
        #number,name,mag,group,redshift,galactic_comp=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
        #del(number,mag,redshift)
        
    #-----------------------------defining fe template fitter ------------------------------------#
        
        
        
        def model(xarr, scale_fe,shift_fe):
                
            try:
		    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr.value+shift_fe, 
									       fe_template.data )
            except:    
		        new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr+shift_fe, 
									       fe_template.data )
            new_template=  new_fe_template 
            return new_fe_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                   2, 
                                                                   parnames=['scale_fe','shift_fe'],#,'scale_balmer'], 
                                                                   parlimited=[(True,True),(True,True)],
                                                                   parlimits=[(slim[0],slim[1]), (vlim[0],vlim[1])],
                                                                   
                                                                    shortvarnames=(r'Af',r'\Delta xf'),#,r'Ab'),
                                                                   centroid_par='shift_fe' )
        fitter_name="rescale"+str(np.random.rand(1)[0])
        modelclass.__name__ = fitter_name 

            

        #-----------------------------defining fe template fitter ------------------------------------#

        
        
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
            
        sp.Registry.add_fitter(fitter_name,modelclass,2,override=True,multisingle='multi')
        sp.specfit(fittype=fitter_name,exclude=exclude,guesses=[1.0,0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
        fe_template.data*=sp.specfit.parinfo.values[0]
        
	try:
		wltemp=fe_template.xarr.value+sp.specfit.parinfo.values[1]
		fe_template.xarr=wltemp
	except:
		fe_template.xarr+=sp.specfit.parinfo.values[1]
        
        
        
        fitmin=np.argmin(np.abs(sp.xarr.value-xmin))
        fitmax=np.argmin(np.abs(sp.xarr.value-xmax))
        print fitmin, fitmax, len(sp.data)
        
        if False:
            #-----------restarting plotter------------#
            sp.plotter.refresh()    
            sp.plotter.reset_limits()
            sp.plotter(figure=1,xmin=xmin,xmax=xmax)
            #-----------restarting plotter------------#
        
            plot_file=plot_path + "fe_fit_" + line_name + "_" + group[spec_number]+"_"+ name[spec_number] + ".png"
            pylab.figure()
            pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
            pylab.xlim(xmin=xmin-300,xmax=xmax)
            pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
            pylab.xlabel(r'$\AA$')
            pylab.plot(sp.xarr,sp.data,'k')
            pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
            pylab.savefig(plot_file)
            
            #pylab.close()
            # -----------------  ploting fit---------------------
        
        sp.specfit.fullsizemodel() 
        sp.data=sp.data-sp.specfit.fullmodel
    return sp.specfit.fullmodel, sp.specfit.parinfo.values





def fe_fitter_shift(sp,line_name,spec_number,xmin,xmax,fe_template, magorder,plot_path="./plots/", do_fit=True):
    if do_fit:
        
	    
        
    #-----------------------------defining fe template fitter ------------------------------------#
        
        
        
        def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                

	                    
            try:
		    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr.value+shift_fe, 
									       fe_template.data )
            except:    
		    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
									       fe_template.xarr+shift_fe, 
									       fe_template.data )

                





            return new_fe_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                       3, 
                                                                       parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                       parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                       parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                         
                                                                       shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                       centroid_par='shift_fe' )
        modelclass.__name__ = "template" 

            

        #-----------------------------defining fe template fitter ------------------------------------#

        
        #-----------restarting plotter------------#
        #sp.plotter.refresh()    
        #sp.plotter.reset_limits()
        #sp.plotter(figure=1,xmin=xmin,xmax=xmax)
        #-----------restarting plotter------------#
        
        
        # -----------------  fitting ---------------------
        
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
        
        sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
        sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
        fe_template.data*=sp.specfit.parinfo.values[0]
        
	try:
		wltemp=fe_template.xarr.value+sp.specfit.parinfo.values[1]
		fe_template.xarr=wltemp
	except:
		fe_template.xarr+=sp.specfit.parinfo.values[1]
        


        
        new_template=fe_template.copy()            
        deltav=100.0/2.355
        sigmatemplate=900.0/2.355
        shift_total=np.arange(-20.0,21.0,1.0)    
            
        templates=np.empty( (len(fe_template.data),100) )
        chi2=np.empty(100)
        amplitude=np.empty(100)
        shift=np.empty(100)
        fwhm=np.empty(100)
        
        templates[:,0]=fe_template.data
        
        chi2[0]=sp.specfit.chi2
        amplitude[0]=sp.specfit.parinfo.values[0]
        shift[0]=sp.specfit.parinfo.values[1]
        fwhm[0]=sigmatemplate*2.355
        
        for i in range(1,100):
            for shift_fe in shift_total:    
                delta_sigma=np.sqrt(2*deltav*sigmatemplate+deltav**2)
                sigma=kms_to_wl(delta_sigma,2800) 
                templates[:,i]=gauss_conv(templates[:,i-1], sigma, order=0, output=None, mode='reflect', cval=0.0)
                fe_template.data=templates[:,i]
                sigmatemplate=sigmatemplate+deltav
            
                def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
                    

                    try:
			    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
										       fe_template.xarr.value+shift_fe, 
										       fe_template.data )
		    except:    
			    new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr, 
										       fe_template.xarr+shift_fe, 
										       fe_template.data )



                    
                    return new_fe_template


        


                modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                           3, 
                                                                           parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                           parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                           parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                           shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                           centroid_par='shift_fe' )
		modelname="template" + str(np.random.rand(1)[0])
                modelclass.__name__ = modelname 


        
            
            
        
                sp.Registry.add_fitter(modelname,modelclass,3,multisingle='multi')
                sp.specfit(fittype=modelname,exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax,debug=False,verbose=False)
            
                sp.specfit.fullsizemodel()  
            
            
            
                chi2[i]=sp.specfit.chi2
                amplitude[i]=sp.specfit.parinfo.values[0]
                shift[i]=sp.specfit.parinfo.values[1]
                fwhm[i]=sigmatemplate*2.355
                
        arg_chi=np.argmin(chi2)
        #arg_chi=0
        chi2_min=chi2[arg_chi]
        amplitude_min=amplitude[arg_chi]
        shift_min=shift[arg_chi]
        fwhm_min=fwhm[arg_chi]
        new_template.data=templates[:,arg_chi]
        new_template.data*=amplitude_min
        new_template.xarr+=shift_min
        out_params=np.empty(3)    
        out_params[0]=fwhm_min
        out_params[1]=shift_min
        out_params[2]=amplitude_min
        
        print np.argsort(chi2)[:10]

        new_template = pyspeckit.interpolation.interp(new_template,sp)
                    
        sp.specfit.residuals=sp.data - new_template.data
        sp.specfit.fullmodel=new_template.data
        
        print "arg_chi,chi,amp,shift,fwhm= ",arg_chi,chi2_min,amplitude_min,shift_min,fwhm_min
        # -----------------  fitting ---------------------
        
    # -----------------  ploting  fit---------------------
        
        sp.specfit.fullsizemodel() 
        fitmin=np.argmin(np.abs(sp.xarr-xmin))
        fitmax=np.argmin(np.abs(sp.xarr-xmax))
        print fitmin, fitmax, len(sp.data)
        
        plot_file=plot_path + "fe_fit_" + line_name   + "_0.png"
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=xmin-300,xmax=xmax)
        pylab.ylabel(r'$10^{'+str(magorder-8)+'}$ erg s$^{-1}$  $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr,sp.data,'k',label="FWHM="+str(fwhm_min))
        pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
        pylab.savefig(plot_file)
        pylab.legend()
        #pylab.close()
        # -----------------  ploting fit---------------------
        
        
       # ----------plotting residuals-------------
        
        pylab.figure()
        pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
        pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
        pylab.ylabel(r'$10^{-'+str(magorder-8)+'}$ erg s$^{-1}$ $\AA^{-1}$')
        pylab.xlabel(r'$\AA$')
        pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
        plot_file=plot_path + line_name +"_res_" + group[spec_number]+"_"+ name[spec_number] + ".png"
        pylab.savefig(plot_file)
        #pylab.close()
        # ----------plotting residuals-------------#
    
        #----------redifining spectrum ------------------#
        sp.data=sp.data - sp.specfit.fullmodel
        #----------redifining spectrum ------------------#
        model=sp.specfit.fullmodel
    else:
        model=0.0*sp.data
    return model,out_params




def convolve_to_fwhm_mg():
    
    sp.data=sp1.data  # DO NOT FORGET TO INCLUDE THIS LINE
    #---------Normalizing Fe Template -------------
    xmin_fe=2200
    xmax_fe=3648
    
    xmin_bb=3640
    xmax_bb=3650
    lambda_0=3675
    index=np.argmin( np.abs( lambda_0 - sp.xarr) )
    flux=np.mean(sp.data[index -30: index + 30])
    template_file=fe_template_UV_file
    fe_template_UV = pyspeckit.Spectrum(template_file)
    fe_template_UV.xarr.unit='angstroms'
    fe_template_UV.xarr.units='angstroms'
    fe_template_UV.xarr.xtype='wavelength'
    #fe_template_UV.crop(1555.0,fe_template_UV.xarr[-1],unit="angstroms")
    lambda_fe=2700
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template_UV.xarr) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template_UV.data*=flux/fe_template_UV.data[index_fe]
    #---------Normalizing Fe Template -------------
    
            
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
            



def load_op_template(sp,mag_order,fe_template_OP_file,galaxy_template_file):

    template_file=fe_template_OP_file
    fe_template = pyspeckit.Spectrum(template_file)
    galaxy_template=pyspeckit.Spectrum(galaxy_template_file)
    xmin_fe_template=fe_template.xarr[0]
    fe_template.xarr.unit==un.Angstrom#'angstroms'
    #fe_template.xarr.unit='angstroms'
    fe_template.xarr.xtype='wavelength'
    #galaxy_template.xarr.unit='angstroms'
    galaxy_template.xarr.unit==un.Angstrom#'angstroms'
    galaxy_template.xarr.xtype='wavelength'
    galaxy_template.data*=1e8
    galaxy_template.data*=10**(-1*mag_order)
    galaxy_template = pyspeckit.interpolation.interp(galaxy_template,sp)
    lambda_fe=4750
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr.value) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template.xarr.value) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template.data*=np.abs(flux/fe_template.data[index_fe])
    return fe_template,galaxy_template

def op_continuum(i,sp,mag_order,OP_limits,w,galactic_comp,plot_objpath): 
    continuous,wlmin_OP,L_model=continuous_substraction( i, sp, mag_order,OP_limits,w)
    cont=sp.copy()
    #cont.data=continuous
    #cont,galactic_cont,host_parameters=quasar_host_continuous(i,sp,cont,np.float(galactic_comp),mag_order,w,plot_path=plot_objpath)
    #continuous=cont.data
    host_parameters=[0,0]
    galactic_cont=0.0*sp.data
    print "host_parameters= ", host_parameters[0],host_parameters[1]
    return continuous, wlmin_OP,L_model,galactic_cont,host_parameters



def load_uv_template(sp,mag_order,fe_template_UV_file):
    

    
    #---------Normalizing Fe Template -------------
    xmin_bb=3640
    xmax_bb=3650
    lambda_0=3640
    index=np.argmin( np.abs( lambda_0 - sp.xarr.value) )
    flux=np.mean(sp.data[index -30: index + 30])
    template_file=fe_template_UV_file

    fe_template_UV = pyspeckit.Spectrum(template_file)
            
    # -----------set up units properly------------------#        
    fe_template_UV.xtype = 'wavelength'
    fe_template_UV.xarr.set_unit(un.Angstrom)
    #sp.xarr.units='angstrom'

    
    #-------------- set up unit properly------------#


    #fe_template_UV.xarr.set_units=un.Angstrom#'angstroms'
    #fe_template_UV.xarr.units='angstroms'
    #fe_template_UV.xarr.xtype='wavelength'
    #fe_template_UV.crop(1555.0,fe_template_UV.xarr[-1],unit="angstroms")
    lambda_fe=2700
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr.value) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template_UV.xarr.value) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template_UV.data*=np.abs(flux/fe_template_UV.data[index_fe])
    #---------Normalizing Fe Template -------------
    return fe_template_UV
            


def uv_continuum(i,sp,mag_order,UV_limits,w,balmer_cont_template,balmer_lines_template,plot_objpath): 
    continuous,wlmin_UV,L_model=continuous_substraction( i, sp, mag_order,UV_limits,w)
    backup=sp.copy()
    sp.data=sp.data - continuous
    balmer_template,balmer_tot, index=balmer_normalization(sp,balmer_cont_template,balmer_lines_template)
    
    sp.data=backup.data
    
    cont=sp.copy()
    cont.data=continuous
    
    cont,balmer_template,scale_balmer=total_continuous_fit(i,sp,cont,balmer_template,mag_order,w,plot_path=plot_objpath)
    continuous=cont.data
    balmer_tot=scale_balmer*balmer_tot
    k=0
    
    #-------extracting balmer continuum (not lines yet)-------#
    sp.data[:index]=sp.data[:index] - balmer_template.data[:index] 
    #-------extracting balmer continuum (not lines yet)-------#
    return  sp,continuous, balmer_tot,balmer_template,L_model,scale_balmer,index,wlmin_UV

