# BLRmodelling
BLR line fitting and property measurement of the near UV and optical most prominent  broad emission lines in type 1 AGN. 

To fit Lya, SiOIV, CIV, CIII and Halpha lines: 


python CIV_CIII_SiOIV_modelling_general.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'
# Please make your spectrum data files in path_to_data/ to have '.txt' extension. 


To fit MgII:

python MgII_modelling_general.py MgII path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'

To fit Hbeta line:

python Hbeta_modelling_general.py Hbeta path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'


To measure emission line properties:

python measuring_quasar_properties_general_errors.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name
# remember to add after 'path_to_data'  the symbol '/'

To do nice plots of the fit line properties (still needs development for MgII and Hbeta):

python plot_best_fits.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name
# remember to add after 'path_to_data'  the symbol '/'




The file constraint_single.cfg is the configuration file, here I define the  emission line complexes and set the limits of the 
parameters associated to the gaussian 
components. Each BLR compoment is determined by two gaussians and each NLR component with one gaussian.

The user should modify the limits as follows:


em_lw=150           # km/s minimum width allowed for gaussian components of the broad emision lines (sigma=FWHM/2.355). 
em_uw=8000          #km/s  maximum width allowed for  gaussian components of broad emision lines
em_uw_hb=8000       #km/s  maximum width allowed for broad emision lines
em_uw_sf=3000       #km/s  maximum width allowed for broad emision lines. This is the maximum with for semiforbidden lines like CIII] 
emn_lw=30           # km/s minimum width allowed for narrow emision lines. This should not be smaller that the spectral bin. deltaLAMBDA*c/LAMBDA 
emn_uw=150          # km/s maximum width allowed for narrow emision lines. It should in principle be smaller or equal to em_lw
xsabs_lw=0          #km/s  minimum width allowed for absortion lines
abs_uw=800          #km/s  maximum width allowed for absortion line
d_center=2000       #km/s red-blue shift allowed to the center of  broad emission line. 
		    #If you have confident redshift estimations, typical limits are 1000km/s. However, this is no typically the case  
d_center_n=1000     #km/s red-blue shift allowed to the center of narrow emission line. 
                    #If you have confident redshift estimations, typical limits are 500km/s. However, this is no typically the case  


# The parameters below are not generally determinant, the code is very wise at setting them properly. I would leave the just like they are defined

em_la=0.05          #Minimun amplitude emission lines compared with the medium flux of the spectra.
em_ua=1.0           #Maximun amplitude emission lines compared with the medium flux of the sp
ab_la=-1.0          #Minimun amplitude emission lines compared with the medium flux of the sp
ab_ua=0.0           #Maximun amplitude emission lines compared with the medium flux of the sp