#!/usr/bin/env python

##################################################################################################
# 	Code designed and implemented by Michael Andreas Imhof
#	2021-02-15
##################################################################################################

#--------------------------------------------------------------------
import os.path

# This script launches the ice flow model script
#----------------------------------------

os.system('./imhof-sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_upsurging_ye2000 -model_choice upsurging ')

#os.system('./imhof-sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_m2_ye2000 -model_choice m2 ')
#os.system('./imhof-sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_muscl_ye2000 -model_choice muscl ')

#os.system('./imhof-sia-upsurging.py -smbm eismint1fm -ib None -ye 20000 -of 100 -exp_name eismint1fm -model_choice upsurging ')
#os.system('./imhof-sia-upsurging.py -smbm eismint1mm -ib None -ye 20000 -of 100 -exp_name eismint1mm -model_choice upsurging ')

