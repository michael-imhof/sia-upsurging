#!/usr/bin/env python

'''
##################################################################################################
#    Code designed and implemented by Michael Andreas Imhof 
#    2021-02-15
##################################################################################################


This script contains an SIA ice flow model employing UPSURGING, described in Appendix C in:

Imhof, Michael Andreas,
"Combined climate-ice flow modelling of the Alpine Ice Field during the Last Glacial Maximum",
ETH Zurich, dissertation no. 27416,
DOI:???,
2021


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

from __future__ import division
import numpy as np
import os.path
import time
from copy import deepcopy


import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages


#################################################################3
#    subroutines for plotting
#################################################################3


#------------------------------------------------------------------
def make_fig_endstate(pdf_name, bedrock, ice, surf_elev):
    # this subroutine creates a figure at the end of the simulation

    hmin=0
    hmax=None

    pdf = PdfPages(pdf_name)

    figa = plt.figure(figsize=(12, 12))

    gs1 = GridSpec(1, 1)
    gs1.update(left=0.0, right=0.9, bottom=0, top=1, wspace=0.0, hspace=0.1)
    ax1 = plt.subplot(gs1[0, 0])

    gs1 = GridSpec(1, 1)
    gs1.update(left=0.91, right=1, bottom=0, top=1, wspace=0.0, hspace=0.1)
    ax1c = plt.subplot(gs1[0, 0])

    # draw data
    ax1.imshow(bedrock, cmap='binary', interpolation='nearest', origin='lower', vmin=0, vmax=4000)
    ice_masked = np.ma.masked_where(ice < 0.001, ice) 
    data_plot = ax1.imshow(ice_masked, cmap='plasma_r', interpolation='nearest', origin='lower',vmin=hmin, vmax=hmax)

    # contour lines
    contour_levels = range(1000, 6000, 200)
    ax1.contour(surf_elev,levels=contour_levels, colors='k', linewidths=[2], alpha=0.4)
    contour_levels = range(1000, 6000, 1000)
    ax1.contour(surf_elev,levels=contour_levels, colors='k', linewidths=[2], alpha=1)


    # colorbar
    cbar = plt.colorbar(data_plot, ax1c)
    cbar.set_label(label='Ice thickness (m)', fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    ax1c.set_aspect(25)


    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)


    # close profile plot figure
    pdf.savefig(figa, bbox_inches='tight') 
    pdf.close()
    plt.close(figa)

    print(' *** end state figure saved: '+pdf_name)


################################################################33
#    subroutines that write data and info to txt files
#################################################################3
def write_simulation_settings_file():
    """writes model settings used for the simulation to 'simulation_setting_info.txt'"""

    output_file_name='simulation_setting_info.txt'
    f=open(output_file_name,'w')
    f.write('Ice physics and dynamics:'+'\n')
    f.write('model_choice = '+model_choice+'\n')
    f.write('g = '+str(g)+' m s-2'+'\n')
    f.write('rho_i = '+str(rho_i)+' kg m-3'+'\n')
    f.write('n = '+str(n)+' []'+'\n')
    f.write('A_ice = '+str(A_ice)+' Pa-n a-1'+'\n')
    f.write('dx = '+str(dx)+' m'+'\n')
    f.write('\n')
    f.write('Massbalance model: '+smbm+'\n')
    if smbm=='elevation':
        f.write('ELA = '+str(ela)+' m'+'\n')
        f.write('mass balance gradient = '+str(mbal_grad)+' a-1'+'\n')
        f.write('maximum accumulation = '+str(max_acc)+' m a-1'+'\n')
    f.write(''+'\n')
    f.close()

    print(' *** text file with simulation settings written')


#------------------------------------------------------------------
def save_simulation_statistics():
    """writes details of computational performance to 'simulation_statistics_info.txt'"""

    output_file_name='simulation_statistics_info.txt'
    f=open(output_file_name,'w')
    f.write('Average model years / h:  '+str((time_now-ys) / ((clock_now - real_time_start ) / 3600))+'\n')
    f.write('Simulation duration   h:  '+str((clock_now - real_time_start) / 3600)+'\n')
    f.write('Number of iterations:     '+str(iter_nr )+'\n')
    f.write('Average time step length: '+str((ye-ys)/iter_nr )+'\n')
    f.close()

    print(' *** text file with simulation statistics written')

#------------------------------------------------------------------
def write_ts(cumul_time_series_file_name,time_now,ice,vol_error_cum):
    """appends time stamp, ice volume and cumulative volume error to time series file"""

    f_ts=open(cumul_time_series_file_name,'a')
    f_ts.write(str(time_now)+'\t'+str( np.sum(np.sum(ice,axis=1),axis=0)*dx*dx )+'\t'+str(vol_error_cum)+'\n')
    f_ts.close()


#------------------------------------------------------------------
def write_data_to_txt(output_file_name, data, nx, ny):
    """writes the ice thickness to a txt file"""

    f=open(output_file_name,'w')
    for iy in range(ny):
        f.write(str(data[0][iy]))
        for ix in range(1,nx):
            f.write('\t'+str(data[ix][iy]))
        f.write('\n')

    f.close()

    print(' *** text file with simulation data written')

################################################################
#    subroutines for the ice flow modelling
#------------------------------------------------------------------
def dirichlet_bc(ice):
    """set ice thickness to 0 at the domain boundary"""
    ice[:, 0] = 0
    ice[:,-1] = 0
    ice[ 0,:] = 0
    ice[-1,:] = 0
    return ice

#------------------------------------------------------------------
def initialize_model(input_bedrock_file, dx):
    """This subroutine initialize arrays for bedrock (bedrock), domain size (nx,ny), resolution (dx), 
    ice thickness (ice), surface elevation (surf_elev), mass balance (smb)"""

    if(input_bedrock_file == 'None'):

        nx = 31
        ny = 31
        dx = 50000
        bedrock= np.zeros((nx,ny))
        ice = np.zeros((nx,ny))
        smb = np.zeros((nx,ny))

    else:

        # load bedrock from txt file
        bedrock = np.genfromtxt(input_bedrock_file, unpack=True)
        nx = bedrock.shape[0]
        ny = bedrock.shape[1]
        ice = np.zeros((nx,ny))
        smb = np.zeros((nx,ny))

        print(' *** Loaded input bedrock from file: '+input_bedrock_file)


    surf_elev = bedrock + ice
    surf_elev[surf_elev<0] = 0




    return bedrock,nx,ny,dx,ice,surf_elev,smb

#------------------------------------------------------------------
def mass_balance(smbm, ice, surf_elev, time_now, ela, mbal_grad, max_acc, nx, ny, dx):
    """This subroutine calculates the surface mass balance"""

    # init empty smb
    smb= np.zeros((nx,ny)) # m ice /a

    if(smbm=='elevation'):
        smb[:,:] = (surf_elev[:,:]-ela)*mbal_grad
        smb[smb[:,:]>max_acc] = max_acc

    elif(smbm=='eismint1fm'):
        smb=smb + 0.3

    elif(smbm=='eismint1mm'):
        s = 1e-2/1000.0
        R_el = 450000
        for ix in range(nx):
            for iy in range(ny):
                smb[ix,iy]=min(0.5, s*(R_el - dx*((15 - ix)**2 + (15 - iy)**2)**0.5))

    return smb

#------------------------------------------------------------------
def phi(r):
    """ calculates phi for the muscl flux limiter method. 
    This subroutine is not my work and stems from: https://github.com/alexjarosch/sia-fluxlim (2021-02-15)"""

    # Koren
    # val_phi = np.fmax(0,np.fmin(np.fmin((2.*r),(2.+r)),2.))
    
    # superbee
    val_phi = np.fmax(0,np.fmin(np.fmin(2.*r,1.),np.fmin(r,2.)))

    return val_phi

#------------------------------------------------------------------
def ice_flow_sia_muscl(H, S, Gamma, n, dx, dy, CFL, dt_uplim, dt_min):
    """This subroutine calculates the new ice thickness with the muscl method. 
    This subroutine stems from: https://github.com/alexjarosch/sia-fluxlim (2021-02-15)
    I adapted variable names such it can be used with the rest of my code"""

    Ny, Nx = np.shape(S)
    
    k = np.arange(0,Ny)
    kp = np.hstack([np.arange(1,Ny),Ny-1])
    kpp = np.hstack([np.arange(2,Ny),Ny-1,Ny-1])
    km = np.hstack([0,np.arange(0,Ny-1)])
    kmm = np.hstack([0,0,np.arange(0,Ny-2)])
    l = np.arange(0,Nx)
    lp = np.hstack([np.arange(1,Nx),Nx-1])
    lpp = np.hstack([np.arange(2,Nx),Nx-1,Nx-1])
    lm = np.hstack([0,np.arange(0,Nx-1)])
    lmm = np.hstack([0,0,np.arange(0,Nx-2)])
    
    # calculate ice thickness
#    H = S-B
    
    ### all the k components
    # MUSCL scheme UP
    r_k_up_m = (H[np.ix_(k,l)]-H[np.ix_(km,l)])/(H[np.ix_(kp,l)]-H[np.ix_(k,l)])
    H_k_up_m = H[np.ix_(k,l)] + 0.5 * phi(r_k_up_m)*(H[np.ix_(kp,l)]-H[np.ix_(k,l)])
    r_k_up_p = (H[np.ix_(kp,l)] - H[np.ix_(k,l)])/(H[np.ix_(kpp,l)] - H[np.ix_(kp,l)])
    H_k_up_p = H[np.ix_(kp,l)] - 0.5 * phi(r_k_up_p)*(H[np.ix_(kpp,l)] - H[np.ix_(kp,l)])
    
    # calculate the slope gradient;
    s_k_grad_up = (((S[np.ix_(k,lp)] - S[np.ix_(k,lm)] + S[np.ix_(kp,lp)] - S[np.ix_(kp,lm)])/(4*dx))**2. + ((S[np.ix_(kp,l)] - S[np.ix_(k,l)])/dy)**2.)**((n-1.)/2.)
    D_k_up_m = Gamma * H_k_up_m**(n+2.) * s_k_grad_up
    D_k_up_p = Gamma * H_k_up_p**(n+2.) * s_k_grad_up
    
    D_k_up_min = np.fmin(D_k_up_m,D_k_up_p)
    D_k_up_max = np.fmax(D_k_up_m,D_k_up_p)
    D_k_up = np.zeros(np.shape(H))
    # include the local slope to identify upstream
    D_k_up[np.logical_and(S[np.ix_(kp,l)]<=S[np.ix_(k,l)],H_k_up_m<=H_k_up_p)] = D_k_up_min[np.logical_and(S[np.ix_(kp,l)]<=S[np.ix_(k,l)],H_k_up_m<=H_k_up_p)]
    D_k_up[np.logical_and(S[np.ix_(kp,l)]<=S[np.ix_(k,l)],H_k_up_m>H_k_up_p)] = D_k_up_max[np.logical_and(S[np.ix_(kp,l)]<=S[np.ix_(k,l)],H_k_up_m>H_k_up_p)]
    D_k_up[np.logical_and(S[np.ix_(kp,l)]>S[np.ix_(k,l)],H_k_up_m<=H_k_up_p)] = D_k_up_max[np.logical_and(S[np.ix_(kp,l)]>S[np.ix_(k,l)],H_k_up_m<=H_k_up_p)]
    D_k_up[np.logical_and(S[np.ix_(kp,l)]>S[np.ix_(k,l)],H_k_up_m>H_k_up_p)] = D_k_up_min[np.logical_and(S[np.ix_(kp,l)]>S[np.ix_(k,l)],H_k_up_m>H_k_up_p)]
    
    # MUSCL scheme DOWN
    r_k_dn_m = (H[np.ix_(km,l)]-H[np.ix_(kmm,l)])/(H[np.ix_(k,l)]-H[np.ix_(km,l)])
    H_k_dn_m = H[np.ix_(km,l)] + 0.5 * phi(r_k_dn_m)*(H[np.ix_(k,l)]-H[np.ix_(km,l)])
    r_k_dn_p = (H[np.ix_(k,l)] - H[np.ix_(km,l)])/(H[np.ix_(kp,l)] - H[np.ix_(k,l)])
    H_k_dn_p = H[np.ix_(k,l)] - 0.5 * phi(r_k_dn_p)*(H[np.ix_(kp,l)] - H[np.ix_(k,l)])
    
    # calculate the slope gradient;
    s_k_grad_dn = (((S[np.ix_(km,lp)] - S[np.ix_(km,lm)] + S[np.ix_(k,lp)] - S[np.ix_(k,lm)])/(4*dx))**2. + ((S[np.ix_(k,l)] - S[np.ix_(km,l)])/dy)**2.)**((n-1.)/2.)
    D_k_dn_m = Gamma * H_k_dn_m**(n+2.) * s_k_grad_dn
    D_k_dn_p = Gamma * H_k_dn_p**(n+2.) * s_k_grad_dn
    
    D_k_dn_min = np.fmin(D_k_dn_m,D_k_dn_p)
    D_k_dn_max = np.fmax(D_k_dn_m,D_k_dn_p)
    D_k_dn = np.zeros(np.shape(H))
    # include the local slope to identify upstream
    D_k_dn[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(km,l)],H_k_dn_m<=H_k_dn_p)] = D_k_dn_min[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(km,l)],H_k_dn_m<=H_k_dn_p)]
    D_k_dn[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(km,l)],H_k_dn_m>H_k_dn_p)] = D_k_dn_max[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(km,l)],H_k_dn_m>H_k_dn_p)]
    D_k_dn[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(km,l)],H_k_dn_m<=H_k_dn_p)] = D_k_dn_max[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(km,l)],H_k_dn_m<=H_k_dn_p)]
    D_k_dn[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(km,l)],H_k_dn_m>H_k_dn_p)] = D_k_dn_min[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(km,l)],H_k_dn_m>H_k_dn_p)]
    
    ### all the l components
    # MUSCL scheme UP
    r_l_up_m = (H[np.ix_(k,l)]-H[np.ix_(k,lm)])/(H[np.ix_(k,lp)]-H[np.ix_(k,l)])
    H_l_up_m = H[np.ix_(k,l)] + 0.5 * phi(r_l_up_m)*(H[np.ix_(k,lp)]-H[np.ix_(k,l)])
    r_l_up_p = (H[np.ix_(k,lp)] - H[np.ix_(k,l)])/(H[np.ix_(k,lpp)] - H[np.ix_(k,lp)])
    H_l_up_p = H[np.ix_(k,lp)] - 0.5 * phi(r_l_up_p)*(H[np.ix_(k,lpp)] - H[np.ix_(k,lp)])
    
    # calculate the slope gradient;
    s_l_grad_up = (((S[np.ix_(kp,l)] - S[np.ix_(km,l)] + S[np.ix_(kp,lp)] - S[np.ix_(km,lp)])/(4*dy))**2. + ((S[np.ix_(k,lp)] - S[np.ix_(k,l)])/dx)**2.)**((n-1.)/2.)
    D_l_up_m = Gamma * H_l_up_m**(n+2.) * s_l_grad_up
    D_l_up_p = Gamma * H_l_up_p**(n+2.) * s_l_grad_up
    
    D_l_up_min = np.fmin(D_l_up_m,D_l_up_p)
    D_l_up_max = np.fmax(D_l_up_m,D_l_up_p)
    D_l_up = np.zeros(np.shape(H))
    # include the local slope to identify upstream
    D_l_up[np.logical_and(S[np.ix_(k,lp)]<=S[np.ix_(k,l)],H_l_up_m<=H_l_up_p)] = D_l_up_min[np.logical_and(S[np.ix_(k,lp)]<=S[np.ix_(k,l)],H_l_up_m<=H_l_up_p)]
    D_l_up[np.logical_and(S[np.ix_(k,lp)]<=S[np.ix_(k,l)],H_l_up_m>H_l_up_p)] = D_l_up_max[np.logical_and(S[np.ix_(k,lp)]<=S[np.ix_(k,l)],H_l_up_m>H_l_up_p)]
    D_l_up[np.logical_and(S[np.ix_(k,lp)]>S[np.ix_(k,l)],H_l_up_m<=H_l_up_p)] = D_l_up_max[np.logical_and(S[np.ix_(k,lp)]>S[np.ix_(k,l)],H_l_up_m<=H_l_up_p)]
    D_l_up[np.logical_and(S[np.ix_(k,lp)]>S[np.ix_(k,l)],H_l_up_m>H_l_up_p)] = D_l_up_min[np.logical_and(S[np.ix_(k,lp)]>S[np.ix_(k,l)],H_l_up_m>H_l_up_p)]
    
    # MUSCL scheme DOWN
    r_l_dn_m = (H[np.ix_(k,lm)]-H[np.ix_(k,lmm)])/(H[np.ix_(k,l)]-H[np.ix_(k,lm)])
    H_l_dn_m = H[np.ix_(k,lm)] + 0.5 * phi(r_l_dn_m)*(H[np.ix_(k,l)]-H[np.ix_(k,lm)])
    r_l_dn_p = (H[np.ix_(k,l)] - H[np.ix_(k,lm)])/(H[np.ix_(k,lp)] - H[np.ix_(k,l)])
    H_l_dn_p = H[np.ix_(k,l)] - 0.5 * phi(r_l_dn_p)*(H[np.ix_(k,lp)] - H[np.ix_(k,l)])
    
    # calculate the slope gradient;
    s_l_grad_dn = (((S[np.ix_(kp,lm)] - S[np.ix_(km,lm)] + S[np.ix_(kp,l)] - S[np.ix_(km,l)])/(4*dy))**2. + ((S[np.ix_(k,l)] - S[np.ix_(k,lm)])/dx)**2.)**((n-1.)/2.)
    D_l_dn_m = Gamma * H_l_dn_m**(n+2.) * s_l_grad_dn
    D_l_dn_p = Gamma * H_l_dn_p**(n+2.) * s_l_grad_dn
    
    D_l_dn_min = np.fmin(D_l_dn_m,D_l_dn_p)
    D_l_dn_max = np.fmax(D_l_dn_m,D_l_dn_p)
    D_l_dn = np.zeros(np.shape(H))
    # include the local slope to identify upstream
    D_l_dn[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(k,lm)],H_l_dn_m<=H_l_dn_p)] = D_l_dn_min[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(k,lm)],H_l_dn_m<=H_l_dn_p)]
    D_l_dn[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(k,lm)],H_l_dn_m>H_l_dn_p)] = D_l_dn_max[np.logical_and(S[np.ix_(k,l)]<=S[np.ix_(k,lm)],H_l_dn_m>H_l_dn_p)]
    D_l_dn[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(k,lm)],H_l_dn_m<=H_l_dn_p)] = D_l_dn_max[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(k,lm)],H_l_dn_m<=H_l_dn_p)]
    D_l_dn[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(k,lm)],H_l_dn_m>H_l_dn_p)] = D_l_dn_min[np.logical_and(S[np.ix_(k,l)]>S[np.ix_(k,lm)],H_l_dn_m>H_l_dn_p)]
    
    # check the cfl condition
    dt_cfl = CFL * min(dx**2.,dy**2.) / max(max(max(abs(D_k_up.flatten())),max(abs(D_k_dn.flatten()))),max(max(abs(D_l_up.flatten())),max(abs(D_l_dn.flatten()))))
    
    # calculate final diffusion term
    div_k = (D_k_up * (S[np.ix_(kp,l)] - S[np.ix_(k,l)])/dy - D_k_dn * (S[np.ix_(k,l)] - S[np.ix_(km,l)])/dy)/dy
    div_l = (D_l_up * (S[np.ix_(k,lp)] - S[np.ix_(k,l)])/dx - D_l_dn * (S[np.ix_(k,l)] - S[np.ix_(k,lm)])/dx)/dx
    
    div_back = div_k + div_l
    
    # calculate new ice thickness  
    dt = min(dt_uplim, max(dt_min, dt_cfl)) 

    H = H + dt*div_back
    ice_new = deepcopy(H)

    # no negative ice criterion 
    ice_new[ice_new<0.0] = 0.0


    return ice_new, dt, np.sum(np.sum(ice_new - H,axis=1),axis=0)*dx*dx

#------------------------------------------------------------------
def ice_flow_sia_m2(ice, surf_elev, omega, n, dx, nx, ny, CFL, dt_uplim, dt_min, use_upsurging):
    """This subroutine calculates the new ice thickness with the m2 or the upsurging method"""

    ice_new     = np.zeros((nx,ny))
    ice_new_raw = np.zeros((nx,ny))

    # initialize stagered gradients 
    h_x = np.zeros((nx+1,ny  ))
    h_y = np.zeros((nx  ,ny+1))

    h_x_c = np.zeros((nx+1,ny  ))
    h_y_c = np.zeros((nx  ,ny+1))

    h_x2 = np.zeros((nx  ,ny+1))
    h_y2 = np.zeros((nx+1,ny  ))

    # initialize stagered diffusivity
    D_x = np.zeros((nx+1,ny  ))
    D_y = np.zeros((nx  ,ny+1))

    # calculate m2 surface gradients
    h_x[1:nx,:] = (surf_elev[1:nx,:] - surf_elev[0:nx-1,:]) / dx
    h_y[:,1:ny] = (surf_elev[:,1:ny] - surf_elev[:,0:ny-1]) / dx

    if (use_upsurging):
        # cap surface differences with the upstream ice thickness (UPSURGING)
        h_x_c[1:nx,:] = np.maximum(np.minimum(h_x[1:nx,:], ice[1:nx,:]/dx), -ice[0:nx-1,:]/dx) 
        h_y_c[:,1:ny] = np.maximum(np.minimum(h_y[:,1:ny], ice[:,1:ny]/dx), -ice[:,0:ny-1]/dx) 
    else:
        h_x_c[1:nx,:] =  h_x[1:nx,:]
        h_y_c[:,1:ny] =  h_y[:,1:ny]

    # calculate surface gradient along grid cell boundaries 
    h_x2[1:nx-1,1:ny  ] = (h_x[1:nx-1,0:ny-1] + h_x[1:nx-1,1:ny] + h_x[2:nx,0:ny-1] + h_x[2:nx,1:ny])/4.0
    h_y2[1:nx  ,1:ny-1] = (h_y[0:nx-1,1:ny-1] + h_y[1:nx,1:ny-1] + h_y[0:nx-1,2:ny] + h_y[1:nx,2:ny])/4.0

    # calculate diffusivity
    D_x[1:nx,:] = omega*((ice[1:nx,:] + ice[0:nx-1,:])/2.0) **(n+2) *( h_x[1:nx,:]**2.0 + (h_y2[1:nx,:])**2.0 )**((n-1)/2)
    D_y[:,1:ny] = omega*((ice[:,1:ny] + ice[:,0:ny-1])/2.0) **(n+2) *( h_y[:,1:ny]**2.0 + (h_x2[:,1:ny])**2.0 )**((n-1)/2)

    # maximum diffusivity
    D_max = max(np.amax(D_x), np.amax(D_y))

    # calculate maximum allowed timestep
    if (D_max == 0.0):
        dt = dt_uplim
    else:
        dt = min(dt_uplim, max(dt_min, CFL*dx**2/D_max )) 

    # update ice thickness
    ice_new_raw[0:nx,0:ny] = ice[0:nx,0:ny] + dt/dx*( D_x[1:nx+1,0:ny] * h_x_c[1:nx+1,0:ny] - D_x[0:nx,0:ny] * h_x_c [0:nx,0:ny]    
                                                    + D_y[0:nx,1:ny+1] * h_y_c[0:nx,1:ny+1] - D_y[0:nx,0:ny] * h_y_c [0:nx,0:ny])

    ice_new = deepcopy(ice_new_raw)


    # no negative ice criterion 
    ice_new[ice_new<0.0] = 0.0

    return ice_new, h_x, h_y, D_x, D_y, dt, np.sum(np.sum(ice_new - ice_new_raw,axis=1),axis=0)*dx*dx


################################################################33


################################################################33
if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(description='Settings:')
    parser.add_argument('-ib', '--input_bedrock', default='None', type=str, help='input bedrock file path')
    parser.add_argument('-ys', '--ys', default=0.0, type=float, help='simulation starting year')
    parser.add_argument('-ye', '--ye', default=20000.0, type=float, help='simulation ending year')
    parser.add_argument('-of', '--output_freq', default=50.0, type=float, help='output rate (a)')
    parser.add_argument('-dx', '--resolution', default=5000.0, type=float, help='lateral resolution (m)')
    parser.add_argument('-smbm', '--smbm', default='eismint1fm', type=str, help='Surface mass balance model', choices=['eismint1fm', 'eismint1mm','elevation'])
    parser.add_argument('-ela', '--ela', default=3000.0, type=float, help='Equilibrium line altitude (m)')
    parser.add_argument('-mbal_grad', '--mbal_grad', default=0.0075, type=float, help='Mass balance gradient ((m/a)/m)')
    parser.add_argument('-max_acc', '--max_acc', default=3.0, type=float, help='Maximum accumulation rate (m/a)')
    parser.add_argument('-exp_name', '--exp_name', default='spam', type=str, help='Experiment name shown on output directory')
    parser.add_argument('-rt_fig', '--run_time_figure', default='True', type=str, help='Print run time figure of the modelled glaciers', choices=['True','False'])
    parser.add_argument('-model_choice', '--model_choice', default='upsurging', type=str, help='Chose SIA integration scheme', choices=['m2', 'upsurging','muscl'])

    args = parser.parse_args()

#    Initialize strings, booleans and variables from input parser
#------------------------------------------------------------------------
    input_bedrock_file = args.input_bedrock
    ys = args.ys
    ye = args.ye
    dx = args.resolution
    smbm = args.smbm
    output_freq = args.output_freq
    ela = args.ela
    mbal_grad = args.mbal_grad
    max_acc = args.max_acc
    experiment_name = args.exp_name
    model_choice = args.model_choice
    run_time_figure = args.run_time_figure


    if run_time_figure == 'True':
        run_time_figure = True
    else:
        run_time_figure = False


    if model_choice == 'm2':
        use_upsurging = False
    elif model_choice == 'upsurging':
        model_choice = 'm2'
        use_upsurging = True


#=============================================================
#    simulation settings
#=============================================================


#=============================================================
#    Constants Ice dynamics
#=============================================================  
    g = 9.81      # m s-2
    rho_i = 910   # kg m-3
    n = 3         # []
    A_ice = 1e-16 # Pa-n a-1

#=============================================================
#    model time_now
#=============================================================
    time_now = ys
    dt_max = 100.0 # years, time steps are no longer than this
    dt_min =  0.0  # years, time steps are at least this long
    CFL = 0.124    # ()

#=============================================================
#    Load input data and initialize bedrock and ice
#=============================================================
    bedrock, nx, ny, dx, ice, surf_elev, smb = initialize_model(input_bedrock_file, dx)

#=============================================================
#    diagnostic variables
#=============================================================
    vol_error_cum = 0     # m-3
    iter_nr = 0

#=============================================================
#    output dir
#=============================================================
    cwd = os.getcwd()
    time_stamp=time.strftime("%Y_%m_%d__%H_%M_%S", time.gmtime())
    output_dir_name = time_stamp+'_'+experiment_name
    os.system('mkdir '+output_dir_name)
    os.chdir(cwd+'/'+output_dir_name)

    write_simulation_settings_file()

#=============================================================
#    pre process
#=============================================================
    omega = 2/(n + 2)*A_ice*(rho_i*g)**n 
    time_next_out = ys + output_freq

#=============================================================
#    Initialize figure and output files
#=============================================================

    # run time figure
    if run_time_figure:
        fig = plt.figure(figsize=(15, 20)) 
        gs1 = GridSpec(1, 2)
        gs1.update(left=0.03, right=0.85, top=0.97, bottom=0.03, wspace=0.05, hspace=0.05)
        ax1  = [plt.subplot(gs1[0, 0])]
        ax1.append( plt.subplot(gs1[0, 1]))
        gs3 = GridSpec(1, 1)
        gs3.update(left=0.86, right=0.91, top=0.97, bottom=0.03, wspace=0.0, hspace=0.0)
        ax1c  = plt.subplot(gs3[0, 0])


    write_data_to_txt('bedrock_'+str(int(round(time_now))).zfill(6)+'.txt', bedrock, nx, ny)
    write_data_to_txt('ice_'+str(int(round(time_now))).zfill(6)+'.txt', ice, nx, ny)

    # set up time series file
    cumul_time_series_file_name='ts_year-a_vol-m3_cumvolerr-m3.txt'
    f_ts=open(cumul_time_series_file_name,'w')
    f_ts.close()
    write_ts(cumul_time_series_file_name, time_now, ice, vol_error_cum)


    #     start the simulation
    #-----------------------------------

    real_time_start = time.time()
    while(round(time_now,5) < ye):
        iter_t0 = time.time()

        # calculate the maximum allowed timestep to hit the time for writing data to the disk.
        dt_uplim = min(dt_max, time_next_out-time_now)

        # calculate and perform ice dynamics
        if model_choice == 'm2': # or upsurging
            ice_dyn_t0 = time.time()
            ice, h_x, h_y, D_x, D_y, dt, vol_error = ice_flow_sia_m2(ice, surf_elev, omega, n, dx, nx, ny, CFL, dt_uplim, dt_min, use_upsurging)
            ice_dyn_t0 =  time.time()-ice_dyn_t0 

        if model_choice == 'muscl':
            ice_dyn_t0 = time.time()
            ice, dt, vol_error = ice_flow_sia_muscl(ice, surf_elev, omega, n, dx, dx, CFL, dt_uplim, dt_min)
            ice_dyn_t0 =  time.time()-ice_dyn_t0 

        vol_error_cum = vol_error_cum + vol_error

        # calculate mass balance
        smb_t0 = time.time() 
        smb = mass_balance(smbm, ice, surf_elev, time_now, ela, mbal_grad, max_acc, nx, ny, dx)
        smb_t0 =  time.time() - smb_t0 

        # apply mass balance
        ice = ice + smb * dt
        ice[ice<0]=0

        # update ice surface elevation
        surfelev_t0 = time.time() 
        surf_elev[:,:] = np.maximum(bedrock[:,:] + ice[:,:], 0)
        surfelev_t0 =  time.time() - surfelev_t0 

        # dirichlet Boundary Condition
        ice = dirichlet_bc(ice)

        # update time
        time_now = time_now + dt
        iter_nr = iter_nr + 1

        # save model results
        if (round(time_next_out, 5) == round(time_now, 5)):

            time_next_out = time_next_out + output_freq

            # write data to scratch
            write_ts(cumul_time_series_file_name, time_now, ice, vol_error_cum)
            write_data_to_txt('ice_'+str(int(round(time_now))).zfill(6)+'.txt', ice, nx, ny)

            # update run time figure
            if run_time_figure:
                # plot map of bedrock and ice thickness
                ax1[0].clear()
                ax1[0].imshow(bedrock, cmap='binary', interpolation='nearest', origin='lower',vmin=0,vmax=4000 )

                ice_masked = np.ma.masked_where(ice < 0.01, ice) 
                data_plot = ax1[0].imshow(ice_masked, cmap='plasma_r', interpolation='nearest', origin='lower', vmin=0, vmax=None)

                contour_levels = range(1000, 6000, 200)
                ax1[0].contour(surf_elev,levels=contour_levels, colors='k', linewidths=[1], alpha=0.4)
                contour_levels = range(1000, 6000, 1000)
                ax1[0].contour(surf_elev,levels=contour_levels, colors='k', linewidths=[1], alpha=1)

                # plot ice thickness color bar
                ax1c.clear()
                cbar = plt.colorbar(data_plot, ax1c)#, extend='max'
                cbar.set_label(label='Ice thickness (m)', fontsize=20)
                cbar.ax.tick_params(labelsize=20)
                ax1c.set_aspect(10)


                # plot map of bedrock and surface mass balance
                ax1[1].clear()
                ax1[1].imshow(bedrock, cmap='binary', interpolation='nearest', origin='lower',vmin=0,vmax=4000)

                smb_masked= np.ma.masked_where(ice < 0.001, smb) 
                ax1[1].imshow(smb_masked, cmap='coolwarm_r', interpolation='nearest', origin='lower', vmin=-0.3, vmax=0.3)

                ax1[0].get_xaxis().set_visible(False)
                ax1[0].get_yaxis().set_visible(False)
                ax1[1].get_xaxis().set_visible(False)
                ax1[1].get_yaxis().set_visible(False)

                ax1[0].set_title('Ice thickness', fontsize=20)
                ax1[1].set_title('Accumulation/Ablation area', fontsize=20)


                plt.draw()
                plt.pause(1e-17)


        clock_now = time.time()

        # print simulation characteristics to screen

        print('')
        print('Model time_now: '+str(time_now))
        print('------------------------')
        print('dt: '+str(dt)+' a')
        print('Volume error: '+str(vol_error)+' m3')
        print('Cum. volume error: '+str(vol_error_cum)+' m3')
        print('Average model years / h:     '+str((time_now-ys)/((clock_now - real_time_start)/3600)))
        print('Diagnostic model years / h:  '+str(dt/((clock_now - iter_t0)/3600)))
        print('Calc. time ice dynamics (s): '+str(ice_dyn_t0))
        print('Calc. time smb (s):          '+str(smb_t0))
        print('Calc. time surf elev (s):    '+str(surfelev_t0))
        print('Calc. time 1 iteration (s):  '+str(clock_now-iter_t0 ))
        print(round(time_next_out,5))
        print(round(time_now,5))





#=============================================================
#    Save end state figure 
#=============================================================

    save_simulation_statistics()

    make_fig_endstate('end_state.pdf', bedrock, ice, surf_elev)


    # leave output directory
    os.chdir(cwd)
    print('\ndone :)\n')

