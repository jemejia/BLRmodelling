# to run python plot_best_fits.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name # remember to after 'path_to_data' to include the symbol '/'
import json
import pylab
import numpy as np
from fitcode.line import *
pylab.rcParams["figure.figsize"]=8,8
if len(sys.argv)<2:
    line_to_fit='CIV'
else:
    line_to_fit=sys.argv[1]


if len(sys.argv)<4:
    data_dir='./data/spectra/' #dir related with the spectra that will be fit
    fn='./data/filenames.txt' #file with the list of name of the spectra located in data_dir
else:
    data_dir=sys.argv[2]
    fn=sys.argv[3]


filenames=np.genfromtxt(fn,unpack=1,dtype='str')

for filename in filenames:
    fp=open(data_dir+ 'fit/'+ filename + '/' + line_to_fit + '_dictionary.json')
    obj_dict=json.load(fp)
    pylab.figure()
    spectrum_file=file=data_dir+filename
    x1,y1=np.genfromtxt(obj_dict['datafile'],unpack=1)
    xmin=x1[0];xmax=x1[-1]
    ymax=1.1*y1.max()
    
    try:
        xobs,yobs,yerrobs=np.genfromtxt(spectrum_file,unpack=1)
    except:
        xobs,yobs=np.genfromtxt(spectrum_file,unpack=1)
    
    mag_order=np.int((1)*np.round(np.log10(np.mean(yobs))))
    yobs = 10**(-1.0*mag_order)*yobs

    xcont,ycont=np.genfromtxt(obj_dict['continuous']['datafile'],unpack=1)
    ycont1=np.interp(xobs,xcont,ycont)
    arg1=np.argmin(np.abs(xobs-xmin))
    arg2=np.argmin(np.abs(xobs-xmax))
    ymin=((yobs-ycont1)[arg1:arg2]).min()
    pylab.plot(xobs,yobs-ycont1,'k',linewidth=3)
    pylab.xlim(xmin,xmax)
    pylab.ylim(ymin,ymax)


    xarr=np.linspace(xmin,xmax,1000)
    x=np.linspace(xmin,xmax,1000)
    y=np.zeros_like(x)
    ybroad=np.zeros_like(x)
    for component in obj_dict['lines'].keys():
        params=obj_dict['lines'][component]['modelpars']
        yarr=gaussian(xarr,params[0],params[1],params[2])
        y=y+yarr
        pylab.plot(xarr,yarr,'k--')
        if component==line_to_fit+'1' or component==line_to_fit+'2':
            ybroad+=yarr
    pylab.plot(x,y,'r',linewidth=1.5)
    pylab.plot(x,ybroad,'b')
    pylab.savefig(data_dir+ 'plots/' + line_to_fit + filename.split('.txt')[0]+'_components.png')
