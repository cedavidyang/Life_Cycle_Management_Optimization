# post processing of simple corrosion results

import os
import numpy as np
import matplotlib as mpl
#mpl.use("pgf")
pgf_with_custom_preamble = {
    "figure.figsize": [3.54, 2.655],
    #"figure.subplot.bottom": 0.14,
    #"figure.subplot.top": 0.93,
    "font.family": "serif", # use serif/main font for text elements
    "font.size": 9, # use font size
    "text.usetex": True,    # use inline math for ticks
    #'text.latex.unicode': True,
    #"pgf.rcfonts": False,   # don't setup fonts from rc parameters
    #"pgf.preamble": [
        #r'\usepackage{fontspec}',         # load additional packages
        #]
}
mpl.rcParams.update(pgf_with_custom_preamble)
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpldatacursor import datacursor

from corrosion.simpleCorrosion import *
from pyre.distributions import Lognormal

def range(smp, percentile1, percentile2):
    return (np.percentile(smp, percentile1), np.percentile(smp, percentile2))

# copy codes from simpleCorrosionLHS
def main(component_type, service_time, icorr_mean, str_yr=None):
    np.random.seed(64)

    if str_yr is None or str_yr < FRP_DESIGN_YR:
        str_yr = 0
    # beam object
    if str_yr == 0:
        rcBeam = generateRcBeamForMC(component_type)
    else:
        frpBeam = generateFrpBeamForMC(component_type)
        frpBeamDesign = generateFrpBeamForDesign(component_type)
        # get initial_design_capacity
        if 'f' in component_type.lower() or 'd' in component_type.lower():
            initial_design_capacity = designFlexStrength(frpBeamDesign, np.array([0]))
        else:
            initial_design_capacity = designShearStrength(frpBeamDesign, np.array([0]))
    # initial variables
    if 'f' in component_type.lower():
        cover_mean = FLEX_COVER_MEAN
        cover_cov = FLEX_COVER_COV
        distance12 = FLEX_DISTANCE12
        region1_no = FLEX_REGION1_NO
        region2_no = FLEX_REGION2_NO
        ds1_mean = FLEX_DS1_MEAN
        ds2_mean = FLEX_DS2_MEAN
        ds1_cov = FLEX_DS1_COV
        ds2_cov = FLEX_DS2_COV
    elif 's' in component_type.lower():
        cover_mean = SHEAR_COVER_MEAN
        cover_cov = SHEAR_COVER_COV
        distance12 = SHEAR_DISTANCE12
        region1_no = SHEAR_REGION1_NO
        region2_no = SHEAR_REGION2_NO
        ds1_mean = SHEAR_DS1_MEAN
        ds2_mean = SHEAR_DS2_MEAN
        ds1_cov = SHEAR_DS1_COV
        ds2_cov = SHEAR_DS2_COV
    elif 'd' in component_type.lower():
        cover_mean = DECK_COVER_MEAN
        cover_cov = DECK_COVER_COV
        distance12 = DECK_DISTANCE12
        region1_no = DECK_REGION1_NO
        region2_no = DECK_REGION2_NO
        ds1_mean = DECK_DS1_MEAN
        ds2_mean = DECK_DS2_MEAN
        ds1_cov = DECK_DS1_COV
        ds2_cov = DECK_DS2_COV

    # percentage array for LHS
    cdf_lhs = np.linspace(0.5/N_LHS_SMP, 1.-0.5/N_LHS_SMP, num=N_LHS_SMP)
    # concrete cover
    concrete_cover = Lognormal('Xc', cover_mean, cover_mean*cover_cov)
    xc1_smp = concrete_cover.rv.ppf(cdf_lhs)
    xc1_smp = np.random.permutation(xc1_smp)
    xc2_smp = xc1_smp + distance12
    # diffusion coefficient
    diffusion_coefficient = Lognormal('Dc', DC_MEAN, DC_MEAN*DC_COV)
    dc_smp = diffusion_coefficient.rv.ppf(cdf_lhs)
    dc_smp = np.random.permutation(dc_smp)
    # surface chloride
    surface_chloride = Lognormal('C0', C0_MEAM, C0_MEAM*C0_COV)
    c0_smp = surface_chloride.rv.ppf(cdf_lhs)
    c0_smp = np.random.permutation(c0_smp)
    # critical chloride
    critical_chloride = Lognormal('Ccr', CCR_MEAN, CCR_MEAN*CCR_COV)
    ccr_smp = critical_chloride.rv.ppf(cdf_lhs)
    ccr_smp = np.random.permutation(ccr_smp)
    # corrosion initiation time
    ti1_smp = xc1_smp**2 / (4.*dc_smp) * (spc.erfinv((c0_smp-ccr_smp)/c0_smp))**(-2)
    ti2_smp = xc2_smp**2 / (4.*dc_smp) * (spc.erfinv((c0_smp-ccr_smp)/c0_smp))**(-2)
    # initial reinforcement
    rebar_region1 = Lognormal('ds1', ds1_mean, ds1_cov)
    ds1_smp = rebar_region1.rv.ppf(cdf_lhs)
    ds1_smp = np.random.permutation(ds1_smp)
    rebar_region2 = Lognormal('ds2', ds2_mean, ds2_cov)
    ds2_smp = rebar_region2.rv.ppf(cdf_lhs)
    ds2_smp = np.random.permutation(ds2_smp)
    resistance_mean = []
    resistance_cov = []
    resistance_smp = []
    # corrosion rate
    corrosion_current = Lognormal('icorr', icorr_mean, icorr_mean*ICORR_COV)
    icorr_smp = corrosion_current.rv.ppf(cdf_lhs)
    icorr_smp = np.random.permutation(icorr_smp)
    corrosion_rate = icorr_smp * 0.0232
    if service_time[0] != 0:
        service_time = np.insert(service_time, 0, 0)
    for indx,t in enumerate(service_time+str_yr):
        # time after corrosion
        t_corrosion1 = t+str_yr - ti1_smp
        t_corrosion1[t_corrosion1<0] = 0
        t_corrosion2 = t+str_yr - ti2_smp
        t_corrosion2[t_corrosion2<0] = 0
        # residual reinforcement
        ds1t_smp = ds1_smp - t_corrosion1*corrosion_rate
        ds1t_smp[ds1t_smp<0] = 0
        ds2t_smp = ds2_smp - t_corrosion2*corrosion_rate
        ds2t_smp[ds2t_smp<0] = 0
        Ast_smp = np.pi/4 * (ds1t_smp**2*region1_no + ds2t_smp**2*region2_no)
        if str_yr == 0:
            rcBeam.Ast = Ast_smp
            rcBeam.Asvt = Ast_smp
            if 'f' in component_type.lower() or 'd' in component_type.lower():
                rsmp = rcBeam.flexCapacity()
            else:
                rsmp = rcBeam.shearCapacity()
            cost = 0
        elif 'f' in component_type.lower():
            #initial_design_capacity = 1.25*M_DCDW_NOM + 1.75*(M_LLIM_NOM+M_DESIGN_LANE)*LDF_FLEX
            frpBeamDesign.ds = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Ast = np.array([np.mean(Ast_smp)])
            frpBeam.ds = ds1t_smp
            frpBeam.Ast = Ast_smp
            if t == str_yr:
                # design FRP strengthening
                design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                    tfrp) - initial_design_capacity
                tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
                if tfrp[0]<0:
                    tfrp = 0.0
                elif tfrp[0]>MAX_LAYER * 0.166:
                    tfrp = MAX_LAYER * 0.166
                else:
                    tfrp = tfrp[0]
                cost = GIRDER_NUM*tfrp*BFRP_NOM*(SPAN_NOM-2*PLATE2SUPPORT_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrp = tfrp * np.ones(N_LHS_SMP)
            rsmp = frpBeam.flexCapacity()
        elif 's' in component_type.lower():
            #initial_design_capacity = 1.25*V_DCDW_NOM + 1.75*(V_LLIM_NOM+V_DESIGN_LANE)*LDF_SHEAR
            frpBeamDesign.dsv = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Asvt = np.array([np.mean(Ast_smp)])
            frpBeam.dsv = ds1t_smp
            frpBeam.Asvt = Ast_smp
            if t == str_yr:
                # design FRP strengthening
                design_frp = lambda tfrpv: designShearStrength(frpBeamDesign,
                    tfrpv) - initial_design_capacity
                tfrpv, infodict, ier, mesg = fsolve(design_frp, 0.001, full_output=True)
                if tfrpv[0]<0:
                    tfrpv = 0.0
                elif tfrpv[0]>MAX_LAYER * 0.166:
                    tfrpv = MAX_LAYER * 0.166
                else:
                    tfrpv = tfrpv[0]
                cost = GIRDER_NUM*tfrpv*WFRPV_NOM*(BEAM_HEIGHT_NOM*2+WEB_WIDTH_NOM)*\
                    (SPAN_NOM/3.*2/SFRPV_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrpv = tfrpv * np.ones(N_LHS_SMP)
            rsmp = frpBeam.shearCapacity()
        elif 'd' in component_type.lower():
            #initial_design_capacity = 1.25*M_DCDW_DECK_NOM + 1.75*(1.2*M_LLIM_DECK_NOM)
            frpBeamDesign.bfrp = frpBeamDesign.b/FRP_DECK_RATIO
            frpBeam.bfrp = frpBeam.b/FRP_DECK_RATIO
            frpBeamDesign.ds = np.array([np.mean(ds1t_smp)])
            frpBeamDesign.Ast = np.array([np.mean(Ast_smp)])
            frpBeam.ds = ds1t_smp
            frpBeam.Ast = Ast_smp
            if t == str_yr:
                # design FRP strengthening
                design_frp = lambda tfrp: designFlexStrength(frpBeamDesign,
                    tfrp) - initial_design_capacity
                tfrp, infodict, ier, mesg = fsolve(design_frp, 0.11, full_output=True)
                if tfrp[0]<0:
                    tfrp = 0.0
                elif tfrp[0]>MAX_LAYER * 0.166:
                    tfrp = MAX_LAYER * 0.166
                else:
                    tfrp = tfrp[0]
                cost = (GIRDER_NUM-1)*tfrp*(SPAN_NOM/FRP_DECK_RATIO)*\
                    (SPAN_DECK_NOM-2*PLATE2SUPPORT_DECK_NOM)*FRP_UNIT_PRICE
                frpBeam.tfrp = tfrp * np.ones(N_LHS_SMP)
            rsmp = frpBeam.flexCapacity()

        rmean = np.mean(rsmp)
        rstd = np.std(rsmp)
        rcov = rstd/rmean
        resistance_mean.append(rmean)
        resistance_cov.append(rcov)
        resistance_smp.append(rsmp)

    return ti1_smp, np.array(resistance_smp), np.array(resistance_mean), np.array(resistance_cov)

if __name__ == '__main__':
    service_time = np.arange(START_AGE+TIME_INTERVAL,END_AGE+TIME_INTERVAL,TIME_INTERVAL)
    if not os.path.isfile('./data/deterioration_bbb.npz'):
        # flexural data
        flex_ti1_smp, resistance_rc, flex_mean, flex_cov = main('flexure', service_time, 1)
        dummy, resistance_frp, dummy, dummy = main('flexure', service_time, 1, 30)
        # shear data
        shear_ti1_smp, dummy, shear_mean, shear_cov = main('shear', service_time, 1)
        # deck data
        deck_ti1_smp, dummy, deck_mean, deck_cov = main('deck', service_time, 1)
    else:
        deterioration = np.load('./data/deterioration_bbb.npz')
        flex_ti1_smp = deterioration['flex_ti1_smp']
        resistance_rc = deterioration['resistance_rc']
        flex_mean = deterioration['flex_mean']
        flex_cov = deterioration['flex_cov']
        resistance_frp = deterioration['resistance_frp']
        shear_ti1_smp = deterioration['shear_ti1_smp']
        shear_mean = deterioration['shear_mean']
        shear_cov = deterioration['shear_cov']
        deck_ti1_smp = deterioration['deck_ti1_smp']
        deck_mean = deterioration['deck_mean']
        deck_cov = deterioration['deck_cov']

    plt.close('all')
    #plt.rc('font', family='serif', size=12)

    # corrosion initiation
    num_bins = 100
    plt.figure()
    plt.hist(flex_ti1_smp, bins=num_bins, # range=range(flex_ti1_smp, 0,100),
            normed=True, color='blue', histtype='step', label='Flexure (girder)')
    plt.hist(shear_ti1_smp, bins=num_bins, # range=range(shear_ti1_smp, 0,100),
            normed=True, color='red', histtype='step', label='Shear (girder)')
    plt.hist(deck_ti1_smp, bins=num_bins, # range=range(deck_ti1_smp, 0,100),
            normed=True, color='green', histtype='step', label='Deck')
    plt.xlim((0, 60))
    plt.xlabel('Time (year)')
    plt.ylabel('PDF')
    annotate_text = u'{label}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    # resistance histogram
    num_bins = 20
    plt.figure()
    resistance_rc_mean = np.mean(resistance_rc[0])
    resistance_rc_std = np.std(resistance_rc[0])
    lognormal_rc = Lognormal('rc', resistance_rc_mean, resistance_rc_std).rv
    resistance_frp_mean = np.mean(resistance_frp[0])
    resistance_frp_std = np.std(resistance_frp[0])
    lognormal_frp = Lognormal('frp', resistance_frp_mean, resistance_frp_std).rv
    plt.hist(resistance_rc[0], bins=num_bins, color='blue', facecolor='none',
            ls='solid', normed=True)
    plt.hist(resistance_frp[0], bins=num_bins, color='red', facecolor='none',
            ls='dashed', normed=True)
    plt.ylim((0,0.002))
    #ax = plt.gca()
    #ax.annotate('initial strength', xy=(1981, 5.e-4), xytext=(2387, 7.29e-4),
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    #ax.annotate('after strengthening', xy=(2079,9.6e-4), xytext=(2387, 1.2e-3),
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    xlims = plt.xlim()
    xticks = np.linspace(xlims[0], xlims[1], 500)
    plt.plot(xticks, lognormal_rc.pdf(xticks), 'b', label='Initial strength')
    plt.plot(xticks, lognormal_frp.pdf(xticks), 'r', label='After strengthening')
    plt.xlabel('Flexural strength of bridge girders (kN-m)')
    plt.ylabel('PDF')
    annotate_text = u'{label}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    # resistance history
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # MC reults
    ax1.plot(service_time, flex_mean[1:]/flex_mean[0], 'b-', label='Flexure (girder)' )
    ax1.plot(service_time, shear_mean[1:]/shear_mean[0], 'r-', label='Shear (girder)' )
    ax1.plot(service_time, deck_mean[1:]/deck_mean[0], 'g-', label='Deck' )
    ax2.plot(service_time, flex_cov[1:], 'b--', label='Flexure (girder)' )
    ax2.plot(service_time, shear_cov[1:], 'r--', label='Shear (girder)' )
    ax2.plot(service_time, deck_cov[1:], 'g--', label='Deck' )
    # figure settings
    ax2.set_xlabel('Time (year)')
    #ax1.set_ylabel('strengh mean (kN-m)')
    #ax2.set_ylabel('strength COV')
    f.text(0.03, 0.5, r'Normalized resistance', ha='center', va='center', rotation='vertical')
    ax1.annotate('Mean values', xy=(5, 0.8), bbox=dict(boxstyle="round", fc="w"))
    ax2.annotate('Coefficient of variation', xy=(5, 0.145), bbox=dict(boxstyle="round", fc="w"))
    ax1.set_ylim((0.75, 1.001))
    #majorFormatter = FormatStrFormatter('%.2e')
    #ax1.yaxis.set_major_formatter(majorFormatter)
    #ax2.yaxis.set_major_formatter(majorFormatter)
    #ax1.annotate('mean (kN-m)', xy=(50.8, 1578), xytext=(19, 1521),
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    #ax2.annotate('COV', xy=(50.4, 0.136), xytext=(20, 0.1369),
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    #ax1.legend(loc='lower left', prop={'size':12})
    #ax2.legend(loc='lower right', prop={'size':12})
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=9,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))


    #plt.show()
