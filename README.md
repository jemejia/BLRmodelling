# BLRmodelling
BLR line fitting and property measurement of the near UV and optical most prominent  broad emission lines in type 1 AGN. 

To fit Lya, SiOIV, CIV, CIII and Halpha lines: 


python CIV_CIII_SiOIV_modelling_general.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'



To fit MgII:

python MgII_modelling_general.py MgII path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'

To fit Hbeta line:

python Hbeta_modelling_general.py Hbeta path_to_data/ file_location_with_list_of_spectra_name 
# remember to add after 'path_to_data'  the symbol '/'


To measure emission line properties:

python measuring_quasar_properties_general_errors.py line_to_fit path_to_data/ file_location_with_list_of_spectra_name
# remember to add after 'path_to_data'  the symbol '/'
