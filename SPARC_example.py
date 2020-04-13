#!/usr/bin/env python3
#
# Version 1.0
#
# This file contains the functions required to genereate Figure 3 of the paper
# "Control of ecological outcomes through deliberate parameter changes in a
# model of the gut microbiome", by Zipeng Wang, Eric Jones, Joshua Mueller, and
# Jean Carlson, published in Physical Review E in 2020. This code implements
# the steady-state reduction (SSR) technique demonstrated in
# steady_state_reduction_example.py (referred to as the SSR code), and applies
# the control framework SPARC to a 11-dimensional generalized Lotka-Volterra
# (gLV) equations fit to C. difficile infection (CDI) gut microbial data. 
#
# In total, this code:
# 1) imports the SSR formulae by Jones and Carlson (Phys. Rev. E 2019)
#    (referred to as the SSR paper) and gLV system parameters fit by Stein et al.
#    (PLOS Comp. Bio. 2013); 
# 2) compresses the high-dimensional parameters into 2D SSR-generated
#    parameters and generate trajectories in both systems starting at the initial
#    condition (0.5,0.5); 
# 3) finds a 2D parameter change that alters the steady-state outcome of this
#    initial condition and plots the resulting new trajectory in the reduced 2D
#    system; and
# 4) generates a high-dimensional parameter modification that corresponds to
#    the 2D parameter change and plots the corresponding trajectory. 
#
# This code is covered under the GNU GPLv3 license.
#
# Send questions to Eric Jones at ewj@physics.ucsb.edu
#
####################################################################################

import steady_state_reduction_example as ssr
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import pickle

####################################################################################
#### HELPER FUNCTIONS
####################################################################################
def integrand(x, t, rho, K):
    """ Integrand for the N-dimensional generalized Lotka-Volterra equations,
    where mu is the growth rates and M is the interaction matrix """
    dxdt = ( np.dot(np.diag(rho), x)
             + np.dot(np.diag(x), np.dot(K, x)) )
    return dxdt

def solve(rho,K,IC,t_end=500,t_num=1001):
    ''' Given a gLV system parametrized by rho (growth rate) and K (interaction
    matrix) and an intitial condition IC, simulates its trajectory in the
    ecological phase space'''
    t = np.linspace(0,t_end, t_num)
    traj = integrate.odeint(integrand,IC,t,args=(rho,K),atol=1e-12)
    return traj

def project_to_2D(traj, ssa, ssb):
    """Projects a high-dimensional trajectory traj into a 2D system,
    defined by the origin and steady states ssa and ssb, and returns a
    2-dimensional trajectory, following Eq. S18 of the supplement of the SSR
    paper. Input must be of list type (even if that means passing "[point]" instead of
    "point"). This function is directl copied from the SSR code"""
    new_traj = []
    for elem in traj:
        uu = np.dot(ssa, ssa); vv = np.dot(ssb, ssb)
        xu = np.dot(elem, ssa); xv = np.dot(elem, ssb)
        uv = np.dot(ssa, ssb)
        new_traj.append([(xu*vv - xv*uv)/(uu*vv - uv**2),
                         (uu*xv - xu*uv)/(uu*vv - uv**2)])
    new_traj = np.array(new_traj)

    return new_traj

def param_change_2D(mu,M,IC):
    '''Given 2D gLV parameters mu and M, generate a change in M_ab that
    changes the steady-state outcome of the initial condition IC (supposing
    that IC tends to the alternative state (0,1) without intervention). This
    function uses the mathematical analysis elaborated in Section II.C.2 '''
    # calculate the nondimensionalized parameters according to Eq. (5)
    nondim_mu_b = mu[1]/mu[0]
    nondim_M_ab = M[0][1]/M[1][1]
    nondim_M_ba = M[1][0]/M[0][0]

    # Since the target state is located at (1,0), we only need to change M_ab
    # The largest change in parameter M_ab that guarantees success 
    # is determined by the bifurcation analysis:
    M_ab_final = M[0][0]*nondim_mu_b

    # the sufficient change of parameter M_ab is then determined by trial and
    # error. For each incremental parameter change, we simulate its trajectory
    # to determine whether it is successful
    eps = 1e-5
    n_steps = 20
    M_ab_step = (M_ab_final - M[0][1])/n_steps
    new_M = M.copy()
    for i in range(n_steps):
        new_M[0][1] += M_ab_step
        trial_traj = solve(mu,new_M,IC)
        if np.linalg.norm(trial_traj[-1]-np.array([1,0]))<eps:
            return new_M

    # If every incremental parameter change fails, it means that the original
    # 2D gLV system is not bistable (i.e., not in the upper-right quadrant of
    # Fig. 2)
    print("2D parameter change failed")

def alpha_gen(ssa,ssb,Change='Mab'):
    ''' given two high-dimensional steady states, and an indicator of which 2D
    interaction matrix entry is being changed, generate the alpha
    coefficients in Eq. (7) according to SSR formula Eq. (3). This function is
    only used in the function parameter_change_11D'''
    ya = ssa/np.linalg.norm(ssa)
    yb = ssb/np.linalg.norm(ssb)
    denom = 1-(sum([ya[k]*yb[k]
               for k in range(len(ya)) ]))**2
    alpha = np.zeros([len(ssa),len(ssa)])
    for i in range(len(ssa)):
        for j in range(len(ssa)):
            if Change == 'Mab':
                numerator = (ya[i]*yb[j]+yb[i]*ya[j])*(ya[i]- yb[i]*sum([ya[k]*yb[k]
                                    for k in range(len(ya))]))
            if Change == 'Mba':
                numerator = (ya[i]*yb[j]+yb[i]*ya[j])*(yb[i]- ya[i]*sum([ya[k]*yb[k]
                                    for k in range(len(ya))]))
            alpha[i][j] = numerator/denom
    return alpha*np.linalg.norm(ssb) # factor comes from rescaling the 2D plane
                                     # accoring to S20-22 of the SSR paper

def parameter_change_11D(K, DeltaM, ssa, ssb, Change='Mab',Entry=1):
    '''Given a high-dimensional interaction matrix K, a two dimensional parameter
    change DeltaM (scalar), and two steady states, generate a modified
    high-dimensional paramter change NewK'''

    alpha = alpha_gen(ssa,ssb,Change)

    # choose the entry K_ij that has the l th largest alpha coefficient
    # (l is given by Entry)
    Flat_alpha=alpha.flatten()
    Flat_alpha.sort()
    target_alpha = Flat_alpha[-Entry]

    # find the coords of the target alpha coefficient
    for i in range(11):
        for j in range(11):
            if alpha[i][j] == target_alpha:
                coord = [i,j]

    # generate NewK
    NewK = K.copy()
    NewK[coord[0]][coord[1]] += DeltaM/target_alpha
    return NewK


def get_lower_bound_of_separatrix(p, s, x_min, x_max, y_min, y_max,
                                  num_sections):
    """ Copied from the SSR code with minor modifications. See
    steady_state_reduction_example.py for a full description """

    delta = (x_max - x_min)/num_sections # spatial variation
    N = num_sections + 1 # number of points needed to make num_sections
    xs = np.linspace(x_min, x_max, num=N)
    num_y_points = int(((y_max - y_min)/delta) + 1)
    ys = np.linspace(y_min, y_max, num=num_y_points)

    ssa = s.ssa
    ssb = s.ssb
    eps = 1e-3
    outcome_dict = {}
    for x in xs:
        print(x)
        for y in ys:
            ic = x*ssa + y*ssb
            y_out = solve(p.rho, p.K, ic, t_end=5000)
            ss = y_out[-1]

            if np.linalg.norm(ss - ssa) < eps:
                outcome_dict[x, y] = 'a'
                #print('a', np.linalg.norm(ss - ssa))
            elif np.linalg.norm(ss - ssb) < eps:
                outcome_dict[x, y] = 'b'
                #print('b', np.linalg.norm(ss - ssb))
            elif x == 0 and y == 0:
                # edge case
                outcome_dict[x, y] = 'a'
            else:
                outcome_dict[x, y] = '?'
                #print('no steady state attained for', x, y)

    lower_bound_list = []
    for x in xs:
        flag = False
        for y in ys[::-1]:
            if flag is False:
                flag = outcome_dict[x, y]
            if flag == outcome_dict[x, y]:
                continue
            if flag != outcome_dict[x, y]:
                lower_bound_list.append([x, y])
                break
                #flag = outcome_dict[x, y]

    return lower_bound_list, delta


####################################################################################
#### FUNCTIONS THAT PRODUCE PLOTS   
####################################################################################
def plot_traj(rho,K,IC,ssa,ssb,ax=None,label=None,color='b'):
    ''' Plots the trajectory in a gLV model parametrized by rho and K with initial condition
    IC '''

    if not ax:
        ax = plt.gca()

    traj = solve(rho,K,IC)

    # If the dimensionality of the model is larger than 2, project the
    # trajectory to 2D and then plot it
    if len(rho)>2:
        traj = project_to_2D(traj,ssa,ssb)

    ax.plot(traj[:,0],traj[:,1],color=color,label= label,lw=2)
    ax.legend() 
    return ax

def plot_2D_separatrix(p, s, ax=None,label=None,color='b'):
    """ Plots the 2D separatrix, as given in Eq 6 of the SSR paper 
    This is directly copied from SSR code with modifications for plotting purposes
    """
    if not ax:
        ax = plt.gca()

    taylor_coeffs = ssr.get_separatrix_taylor_coeffs(s, order=5, dir_choice=1)
    # (u, v) is the coexistent fixed point (x_a^*, x_b^*) of the 2D system
    u, v = s.get_11_ss()
    xs = np.linspace(0, 1, 100)
    ys = [sum( [float(taylor_coeffs[j])*(xx-u)**j/math.factorial(j)
                for j in range(len(taylor_coeffs))] ) for xx in xs]

    ax.plot(xs, ys, '--', color=color, label=label)
    ax.legend(fontsize=12)
    #plt.savefig('SSR_demo_2.pdf')
    #print('... saved figure to SSR_demo_2.pdf')

    return ax

def plot_simple_ND_separatrix(p, s, ax=None, label=None, color='b',
                              sep_filename ='11D_separatrix_1e-2.data'):
    ''' plot the high-dimensional separatrix by sampling points in the bistable
    region'''
    load_data=False
    if load_data:
        with open(sep_filename, 'rb') as f:
            separatrix_lower_bound = pickle.load(f)
    else:
        lower_bound_list, delta = (
            get_lower_bound_of_separatrix(p, s, x_min=0, x_max=1, y_min=0, y_max=1,
                                          num_sections = num_sections) )

        separatrix_lower_bound = np.array(lower_bound_list)
        with open(sep_filename, 'wb') as f:
            pickle.dump(separatrix_lower_bound, f)
        print('(reduce computation time by setting load_data = True in function call for plot_simple_ND_separatrix (line 252)')

    ax.plot(separatrix_lower_bound[:,0], separatrix_lower_bound[:,1],
            '--', color=color, label=label)
    if label:
        ax.legend()
    return ax


def plot_all_panels(p,s,IC):
    '''Given a high-diemnsional gLV model, its SSR-generated 2D model and an
    initial condition of interest, apply the SPARC method and generate Fig. 3
    of the main text. ax[0,0] (upper left) shows the trajectory and the
    separatrix of the original model; ax[1,0] (lower left) shows the
    SSR-generated 2D model; ax[1,1] (lower right) shows the 2D model after
    modification; and ax[0,1] (upper right) shows the high-dimensional model
    after parameter modification.'''
    # Figure configuration
    fig, ax = plt.subplots(2,2)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    fig.set_figheight(10)
    fig.set_figwidth(12)

    font = 16
    xmin,xmax = -0.0,1.05
    ymin,ymax = -0.0,1.05
    for i in range(2):
        for j in range(2):
            ax[i,j].axis([xmin, xmax, ymin, ymax])
            ax[i,j].scatter([1,0,.5],[0,1,.5],facecolor =
                        ['g','r','b'],zorder = 5,clip_on=False)
            #ax[0,1].axis([0,4,0,4])
            if i == 0:
                ax[i,j].set_xlabel(r'$y_a$',fontsize = font)
                ax[i,j].set_ylabel(r'$y_b$',fontsize = font)
            else:
                ax[i,j].set_xlabel(r'$x_a$',fontsize = font)
                ax[i,j].set_ylabel(r'$x_b$',fontsize = font)

            ax[i,j].set_xticks([0,0.5,1])
            ax[i,j].set_yticks([0,0.5,1])

            for side in ['top','bottom','left','right']:
                ax[i,j].spines[side].set_linewidth(2)

    # unpack essential parameters from container class p and s
    rho, K = p.rho, p.K; ssa, ssb = s.ssa, s.ssb; mu, M = s.mu, s.M;

    # (top-left) plot the original high-dimensional trajectory
    IC_highD = IC[0]*ssa + IC[1]*ssb
    plot_traj(rho,K,IC_highD,ssa,ssb,ax[0,0],label='old trajectory')

    # (top-left) plot the original separatrix
    ssr.plot_ND_separatrix(p,s,ax[0,0], label='old separatrix',
                           sep_filename='11D_original_sep.data', delta=0.05)

    # (bottom-left) plot the 2D trajectory in the SSR-generated 2D model
    plot_traj(mu,M,IC,ssa,ssb,ax[1,0],label='old trajectory')

    # (bottom-left) plot the 2D separatrix in this 2D model
    plot_2D_separatrix(p, s, ax[1,0], 'old separatrix')

    # generate the new 2D interaction matrix new_M by trial and error
    new_M = param_change_2D(mu,M,IC)

    # (bottom-right) plot new trajectory of the 2D model after the parameter change
    color_list = [(0,0,1),(0.33,0.16,0.66),(0.66,0.32,0.33),(1,0.55,0)]
    plot_traj(mu,new_M,IC,ssa,ssb,ax[1,1],label='new trajectory',
              color=color_list[-1])

    # (bottom-right) plot the new 2D separatrix four incremental parameter steps
    delta_M_ab = new_M[0][1]-M[0][1]

    steps = 4
    for i in range(steps):
        temp_M = M.copy()
        temp_M[0][1] += (i+1)*delta_M_ab/steps
        temp_s = ssr.ssrParams(p, ssa, ssb)
        temp_s.M = temp_M
        if i==0: i_label = 'old separatrix'
        if i==1 or i==2: i_label = 'incremented separatrix'
        if i==3: i_label = 'new separatrix'

        plot_2D_separatrix(p, temp_s, ax[1,1],label=i_label,
                           color=color_list[i])

    # map this 2D parameter change back to 11D model
    new_K = parameter_change_11D(K,delta_M_ab,ssa,ssb)
    plot_traj(rho,new_K,IC_highD,ssa,ssb,ax[0,1],label='new trajectory',
              color=color_list[-1])

    # (upper-right) plot stepwise separatrices in the upper right panel
    for i in range(steps):
        if i==0: i_label = 'old separatrix'
        if i==1 or i==2: i_label = 'incremented separatrix'
        if i==3: i_label = 'new separatrix'

        temp_K = parameter_change_11D(K, (i+1)*delta_M_ab/steps, ssa, ssb)
        temp_p = ssr.Params()
        temp_p.K = temp_K

        ssr.plot_ND_separatrix(temp_p,s,ax[0,1],sep_filename=f'11D_step{i}_sep.data',
                           color=color_list[i], delta=0.05, label=i_label)


    plt.savefig('SPARC_demo.pdf')
    print('... saved figure to SPARC_demo.pdf')


####################################################################################
#### MAIN FUNCTION   
####################################################################################
if __name__ == '__main__':
    # get Stein Parameters
    p = ssr.Params()

    # generate steady states (labeled A-E) of the Stein model
    stein_ss = ssr.get_all_stein_steady_states(p)

    # designate state E as the target state and state C as the 
    ssa = stein_ss['E']; ssb = stein_ss['C']
    # generate SSR params (mu, M) using steady states E and C
    s = ssr.ssrParams(p,ssa,ssb)
    #mu, M = ssr.get_SSR_params(p, ssa, ssb)
    # Set the initial condition at (0.5,0.5)
    IC = [0.5,0.5]

    # plot Figure 3
    plot_all_panels(p,s,IC)





