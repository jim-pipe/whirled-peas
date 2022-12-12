#         NO CLINICAL USE.  THE SOFTWARE IS NOT INTENDED FOR COMMERCIAL
# PURPOSES AND SHOULD BE USED ONLY FOR NON-COMMERCIAL RESEARCH PURPOSES.  THE
# SOFTWARE MAY NOT IN ANY EVENT BE USED FOR ANY CLINICAL OR DIAGNOSTIC
# PURPOSES.  YOU ACKNOWLEDGE AND AGREE THAT THE SOFTWARE IS NOT INTENDED FOR
# USE IN ANY HIGH RISK OR STRICT LIABILITY ACTIVITY, INCLUDING BUT NOT LIMITED
# TO LIFE SUPPORT OR EMERGENCY MEDICAL OPERATIONS OR USES.  LICENSOR MAKES NO
# WARRANTY AND HAS NOR LIABILITY ARISING FROM ANY USE OF THE SOFTWARE IN ANY
# HIGH RISK OR STRICT LIABILITY ACTIVITIES.
#

# Author: Jim Pipe
# Date: 2022Nov

######################################################
# GPI node to generate WHIRL spiral waveforms
#
# Don't use GPI?  You have two options:
#
# 1.  Go to gpilab.com and become an awesome GPI user
#
# 2.  start below in the section that starts with:
#     def compute(self)
#     and you can use nearly all the python code as is:
#     You will just have to convert the calls that get and pass data
#     This is mostly any line containing the phrase self.xxx
#
# I plan to submit an MRM paper soon, and if accepted I will update this code
# with references to the relevant equations
# No pressure on the reviewers to accept it, though.
######################################################

import gpi

class ExternalNode(gpi.NodeAPI):

    """Module to generate WHIRL (involute of a circle) trajectory and
    associated gradient waveforms using WHIRLED PEAS algorithm.  Also 
    generates associated SDC and time map for image construction.

    OUTPUTS:
    crds_out - output coordinates: the last dimension is 2 (kx/ky).
    grd_out - gradient waveforms used to produce crds_out.
    sdc_out - sampling density weights corresponding to crd array
    time_out - time map of k-space, for use in generating blurring kernels

    WIDGETS: self-explanatory
    Info Box gives the time for 4 WHIRL segments and respective constraints
    ** if constraint is less than what is requested, this was necessary to prevent a
       negative timing - it didn't slow anything down, it's just the max that was needed **
    """

    def execType(self):
        return gpi.GPI_PROCESS

    def initUI(self):

        # Widgets
        self.addWidget('PushButton', 'compute', toggle=True)
        self.addWidget('TextBox', 'Info:')

        self.addWidget('DoubleSpinBox', 'FOV (cm)',
                       val=24.0, min=0.1, decimals=6)
        self.addWidget('DoubleSpinBox', 'Res (mm)',
                       val=0.8, min=0.1, singlestep=0.1, decimals=5)
        self.addWidget('SpinBox', '# of Spiral Arms', val=16, min=1)

        self.addWidget('DoubleSpinBox', 'MaxSlw (mT/m/ms)',
                       val=150.0, min=0.01, decimals=6)
        self.addWidget('DoubleSpinBox', 'MaxGrd (mT/m)',
                       val=40.0, min=0.01, decimals=6)
        self.addWidget('DoubleSpinBox', 'Max G Freq (kHz)', val=1.0, min=0.1)

        self.addWidget('DoubleSpinBox', 'AD dwell time (us)',
                       val=1.0, min=0.1, decimals=6)
        self.addWidget('DoubleSpinBox', 'Grad dwell time (us)',
                       val=1.0, min=0.1, decimals=6)

        self.addWidget('DoubleSpinBox', 'Gam (kHz/mT)',
                       val=42.577, min=0.01, decimals=6, visible = False) # hide this until we need it

        # IO Ports
        self.addOutPort('crds_out', 'NPYarray')
        self.addOutPort('grd_out', 'NPYarray')
        self.addOutPort('sdc_out', 'NPYarray')
        self.addOutPort('time_out', 'NPYarray')

    def validate(self):

        import numpy as np

    def compute(self):

        import numpy as np
        import math

        # convert units to ms, kHz, m, mT
        # pay attention to radial frequencies (omegas) which include a 2Pi factor
        fov = 0.01 * self.getVal('FOV (cm)')
        res = 0.001 * self.getVal('Res (mm)')
        narms = self.getVal('# of Spiral Arms')

        m_slew = self.getVal('MaxSlw (mT/m/ms)')
        m_grad = self.getVal('MaxGrd (mT/m)')
        m_omega = 2.*np.pi*self.getVal('Max G Freq (kHz)')

        dwell = 0.001 * self.getVal('AD dwell time (us)')
        grast = 0.001 * self.getVal('Grad dwell time (us)')
        gamma = self.getVal('Gam (kHz/mT)')

        delta = float(narms)/(2.*np.pi*fov)

        ######################################################################
        # We define a spiral circle to have a diameter (2/sqrt(Pi))*(1/res) in k-space to match a square
        # sampling region of width (1/res), to keep resolution comparable - we call this True Resolution
        #
        # Van Gelderen P. Comparing true resolution in square versus circular k-space sampling.
        # Proceedings of the 6th Annual Meeting of ISMRM, Sydney, Australia, 1998. p 424.
        #
        # We map this to a k-space grid that is 25% larger than (FOV/RES), corresponding to 80% requested resolution
        # This all works out nicely :-) but look down 4 lines if you don't want to do this
        ######################################################################
        trures = math.sqrt(np.pi)/2.
        gridres = 0.8

        # If you want to cheat, you can ignore this by uncommenting the next 2 lines
        # trures = 1.
        # gridres = 1.

        krad_max = 0.5/(trures*res)
        mtx = int(fov/(gridres*res))
        
        ######################################
        # Find compatible constraints, so each segment does not have "negative" duration
        ######################################
        omega1 = math.sqrt(2.*gamma*m_slew/(3.*delta))
        omega2 = 2.*gamma*m_grad/(3.*delta)
        if m_omega > 0:
          omega_max = min([m_omega, omega1, omega2])
        else:
          omega_max = min([omega1, omega2])

        slew1 = m_grad*omega_max
        slew2 = math.sqrt((omega_max**4.)*(krad_max*krad_max - delta*delta)/(gamma*gamma))
        slew_max = min([m_slew,slew1, slew2])

        grad1 = ((slew_max*slew_max)*(krad_max*krad_max - delta*delta)/(gamma*gamma))**0.25
        grad_max = min([m_grad,grad1])

        ######################################
        # Find timings
        # segments start at tx0, end at tx1, with total time t_X
        ######################################
        # Arc
        ta0 = 0
        ta1 = (5*np.pi+1)/(6.*omega_max)
        t_arc = ta1-ta0

        # Omega Constrained
        tw0 = 1./omega_max
        tw1 = (gamma*slew_max)/(delta*omega_max**3.)
        t_omega = tw1-tw0

        # Slew Constrained
        ts0 = (2.*gamma*slew_max)/(3.*delta*omega_max**3.)
        ts1 = (2.*gamma*grad_max**3.)/(3.*delta*slew_max*slew_max)
        t_slew = ts1-ts0

        # Gradient Constrained
        tg0 = (gamma*grad_max**3.)/(2.*delta*slew_max*slew_max)
        tg1 = (krad_max*krad_max - delta*delta)/(2.*gamma*delta*grad_max)
        t_grad = tg1-tg0

        # gradient rampdown is a Hanning Window
        # This may help a little with spiral-in (??)
        t_ramp = np.pi*grad_max/m_slew

        tau_total = t_arc + t_omega + t_slew + t_grad
        tgd_total = tau_total + t_ramp

        ##########################
        # A few things to report back to GPI Info Box
        ##########################
        arc_info = "Arc:  "+format(t_arc,'.3f')+" ms\n"
        frq_info = "Freq: "+format(t_omega,'.3f')+" ms at "+format(omega_max/(2.*np.pi),'.2f')+" kHz \n"
        slw_info = "Slew: "+format(t_slew,'.3f')+" ms at " +format(slew_max,'.2f')+" mT/m/ms \n"
        grd_info = "Grad: "+format(t_grad,'.3f')+" ms at " +format(grad_max,'.2f')+" mT/m \n"
        tau_info = "\nTau: "+format(tau_total,'.3f')+" ms \n"
        tgd_info = "Tgrad: "+format(tgd_total,'.3f')+" ms \n"

        ##########################
        # Compute waveforms
        ##########################
        if self.getVal('compute'):

          gpts = int(tau_total//grast)
          kpts = int(tau_total//dwell)
          rpts = int(t_ramp//grast)

          # k-space points
          krad_out = np.zeros((narms,kpts,2))
          # gradient waveforms
          grad_out = np.zeros((narms,gpts+rpts,2))
          # sampling density calculation
          sdc_out = np.zeros((narms,kpts))
          # Time maps of data acquisition in k-space
          # We use these for a forward model of off-resonance phase accumulation
          # for deblurring, but you may not care about this if you do something else
          time_out = np.zeros((mtx, mtx))

          ##########################
          # Define some constants
          ##########################
          arc_ta = np.pi/(3.*omega_max)
          arc_tb = (1+2.*np.pi)/(6.*omega_max)

          cga = delta*omega_max/(3.*gamma)
          cgw = (delta*omega_max*omega_max/gamma)
          cgs = (3.*delta*slew_max*slew_max/(2.*gamma))**(1./3.)
          cgg = grad_max

          cta = omega_max/3.
          ctw = omega_max
          cts = (9.*gamma*slew_max/(4.*delta))**(1./3.)
          ctg = math.sqrt(2.*gamma*grad_max/delta)

          csa = cga/grad_max
          csw = cgw/grad_max
          css = cgs/grad_max
          csg = cgg/grad_max

          #######################
          # Compute GRADIENT
          #######################
          for i in range(gpts):
            t = float(i)*grast

            #-----
            # ARC |
            #-----
            if t < ta1:
              if t<arc_ta:
                gmag=cga*(1-math.cos(3.*omega_max*t))
                theta = cta*(t-(math.sin(3.*omega_max*t)/(3.*omega_max)))
                theta = theta + 1 - (0.5*np.pi)
              elif t<(arc_ta+arc_tb):
                tt = t - arc_ta
                gmag=2.*cga
                theta = cta*(arc_ta+2.*tt)
                theta = theta + 1 - (0.5*np.pi)
              else:
                tt = t - arc_ta - arc_tb
                gmag=cga*(3-math.cos(3.*omega_max*tt))
                theta = cta*(arc_ta+(2.*arc_tb) + 3.*tt - (math.sin(3.*omega_max*tt)/(3.*omega_max)))
                theta = theta + 1 - (0.5*np.pi)
            else:
              t = t-ta1+tw0

              #-----
              # FRQ |
              #-----
              if t < tw1:
                gmag = cgw*t
                theta = ctw*t
              else:
                t = t-tw1+ts0

                #------
                # SLEW |
                #------
                if t < ts1:
                  gmag = cgs*(t**(1./3.))
                  theta = cts*(t**(2./3.))
                else:
                  t = t-ts1+tg0

                  #------
                  # GRAD |
                  #------
                  if t < tg1:
                    gmag = cgg
                    theta = ctg*math.sqrt(t)
              
            grad_out[0,i,0] = gmag*math.cos(theta)
            grad_out[0,i,1] = gmag*math.sin(theta)

          #############################
          # Compute GRADIENT RAMPDOWN #
          #############################

          for i in range(gpts, gpts+rpts):
            t = float(i)*grast-ta1+tw0-tw1+ts0-ts1+tg0
            theta = ctg*math.sqrt(t)
            t = float(i-gpts)*grast
            gmag = cgg*0.5*(1.+math.cos(np.pi*t/t_ramp))
            grad_out[0,i,0] = gmag*math.cos(theta)
            grad_out[0,i,1] = gmag*math.sin(theta)

          #######################
          # Compute KSPACE and SDC 
          #######################
          for i in range(kpts):
            t = float(i)*dwell
            gmag = cgg
            theta = ctg*math.sqrt(t)

            #-----
            # ARC |
            #-----
            if t < ta1:
              if t<arc_ta:
                sdc=csa*(1-math.cos(3.*omega_max*t))
                theta = cta*(t-(math.sin(3.*omega_max*t)/(3.*omega_max)))
                krad = delta*math.sqrt(2.*(1-math.cos(theta)))
                phi = math.atan2(1-math.cos(theta),math.sin(theta))
                phi = phi+1-0.5*np.pi
                arm2arm = math.sin(theta)
                sdc=csa*(1-math.cos(3.*omega_max*t))*arm2arm
              elif t<(arc_ta+arc_tb):
                tt = t - arc_ta
                theta = cta*(arc_ta+2.*tt)
                krad = delta*math.sqrt(2.*(1-math.cos(theta)))
                phi = math.atan2(1-math.cos(theta),math.sin(theta))
                phi = phi+1-0.5*np.pi
                arm2arm = math.sin(theta)
                sdc=2.*csa*arm2arm
              else:
                tt = t - arc_ta - arc_tb
                theta = cta*(arc_ta+(2.*arc_tb) + 3.*tt - (math.sin(3.*omega_max*tt)/(3.*omega_max)))
                krad = delta*math.sqrt(2.*(1-math.cos(theta)))
                phi = math.atan2(1-math.cos(theta),math.sin(theta))
                phi = phi+1-0.5*np.pi
                arm2arm = math.sin(theta)
                sdc=csa*(3-math.cos(3.*omega_max*tt))*arm2arm
            else:
              t = t-ta1+tw0

              #-----
              # FRQ |
              #-----
              if t < tw1:
                sdc = csw*t
                theta = ctw*t
                krad = delta*math.sqrt(theta*theta+1)
                phi = theta - np.arccos(delta/krad)
              else:
                t = t-tw1+ts0

                #------
                # SLEW |
                #------
                if t < ts1:
                  sdc = css*(t**(1./3.))
                  theta = cts*(t**(2./3.))
                  krad = delta*math.sqrt(theta*theta+1)
                  phi = theta - np.arccos(delta/krad)
                else:
                  t = t-ts1+tg0

                  #-----
                  # GRAD |
                  #-----
                  if t < tg1:
                    sdc = csg
                    theta = ctg*math.sqrt(t)
                    krad = delta*math.sqrt(theta*theta+1)
                    phi = theta - np.arccos(delta/krad)
              
            krad_out[0,i,0] = krad*math.cos(phi)
            krad_out[0,i,1] = krad*math.sin(phi)
            sdc_out[0,i] = sdc

          ########################################
          # calculate SNR factor numerically :-/ #
          # This is the SNR loss by not weighting all data uniformly
          ########################################

          snr_factor = np.mean(sdc_out[0,:])/np.sqrt(np.mean(sdc_out[0,:]*sdc_out[0,:]))
          snr_info = "\nSNR Factor: "+format(snr_factor,'.3f')+"\n"

          #########################
          # Normalize K-space, then
          # Rotate Gradients and K-space to other arms
          #########################
          krad_out[0,:,:] *= gridres*res # output k-space between -0.5 and +0.5 (our standard, but you scale as desired)
          for i in range(1,narms):
            beta = 2.*np.pi*float(i)/float(narms)
            grad_out[i,:,0] = math.cos(beta)*grad_out[0,:,0] - math.sin(beta)*grad_out[0,:,1]
            grad_out[i,:,1] = math.cos(beta)*grad_out[0,:,1] + math.sin(beta)*grad_out[0,:,0]
            krad_out[i,:,0] = math.cos(beta)*krad_out[0,:,0] - math.sin(beta)*krad_out[0,:,1]
            krad_out[i,:,1] = math.cos(beta)*krad_out[0,:,1] + math.sin(beta)*krad_out[0,:,0]
            sdc_out[i,:]    = sdc_out[0,:]

          ######################################
          # Create a Time map for blur kernels #
          # We use these for deblurring        #
          ######################################

          Rc = 0 # just define it for now
          mtx2 = mtx//2
          rc_max = np.sqrt(krad_max*krad_max - delta*delta)

          # Temporal offsets to make smooth transitions
          tbase_frq = ((1+5*np.pi)/(6*omega_max)) - (1./omega_max)
          tbase_slw = tbase_frq + ((gamma*slew_max)/(3.*delta*omega_max*omega_max*omega_max))
          tbase_grd = tbase_slw + ((gamma*grad_max*grad_max*grad_max)/(6.*delta*slew_max*slew_max))

          # Make sure we catch outer edges of time_out
          time_out[:,:] = tau_total

          for i in range(mtx2-1):
            for j in range(mtx2-1):
              rad = np.sqrt(float(i*i+j*j))/fov
              if (rad > delta):
                Rc = np.sqrt(rad*rad-delta*delta)

              #-----
              # Arc |
              #-----
              if (rad < np.sqrt(2.)*delta):
                tacq = (np.pi/(6.*omega_max))+((4*np.pi+1.)*rad/(6.*np.sqrt(2.)*delta*omega_max))

              #-----
              # FRQ |
              #-----
              elif (Rc < (gamma*slew_max)/(omega_max*omega_max)):
                tacq = tbase_frq + (Rc/(delta*omega_max))

              #------
              # SLEW |
              #------
              elif (Rc < (gamma*grad_max*grad_max)/(slew_max)):
                tacq = tbase_slw + (2./(3.*delta))*np.sqrt(Rc*Rc*Rc/(gamma*slew_max))

              #------
              # GRAD |
              #------
              elif (Rc < rc_max):
                tacq = tbase_grd + (Rc*Rc/(2.*gamma*delta*grad_max))

              else:
                tacq = tau_total

              time_out[mtx2+i, mtx2+j] = tacq
              time_out[mtx2+i, mtx2-j] = tacq
              time_out[mtx2-i, mtx2+j] = tacq
              time_out[mtx2-i, mtx2-j] = tacq

          ###########################
          # Report back to info box #
          ###########################
          info = arc_info+frq_info+slw_info+grd_info+tau_info+tgd_info+snr_info
          self.setAttr('Info:',val=info)

          self.setData('crds_out', krad_out)
          self.setData('grd_out', grad_out)
          self.setData('sdc_out', sdc_out)
          self.setData('time_out', time_out)

        return(0)
