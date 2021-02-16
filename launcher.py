#!/usr/bin/env python

'''
##################################################################################################
# 	Code designed and implemented by Michael Andreas Imhof 
#	2021-02-15
##################################################################################################


This script contains exemples for how to launch the main script "sia-upsurging"

Copyright (C) 2021 Michael Andreas Imhof


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


Michael Andreas Imhof was funded by the Swiss National Science Foundation (project 200021-162444).

'''

#--------------------------------------------------------------------
import os.path

# Exemple commands
#----------------------------------------

os.system('./sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_upsurging_ye2000 -model_choice upsurging ')

#os.system('./sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_m2_ye2000 -model_choice m2 ')
#os.system('./sia-upsurging.py -smbm elevation -ib 1km-rhone_valley.nc -dx 1000 -ye 200 -of 10 -ela 2100 -mbal_grad 0.0075 -max_acc 2.00 -exp_name 1km_rhonevalley_muscl_ye2000 -model_choice muscl ')

#os.system('./sia-upsurging.py -smbm eismint1fm -ib None -ye 20000 -of 100 -exp_name eismint1fm -model_choice upsurging ')
#os.system('./sia-upsurging.py -smbm eismint1mm -ib None -ye 20000 -of 100 -exp_name eismint1mm -model_choice upsurging ')

