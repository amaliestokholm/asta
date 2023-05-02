import numpy as np
import uncertainties
import uncertainties.umath
import basta.constants as cs
from basta import seismic_utils as su
from basta.Asfgrid import asfgrid
from . import solar


def calc_corrected_dnufit(osckey, osc, numax=933):
    _, modl0 = su.get_givenl(l=0, osc=osc, osckey=osckey)
    xfitdnu = np.arange(0, len(modl0[0]))
    yfitdnu = modl0[0]
    wfitdnu = np.exp(-1.0 * np.power(yfitdnu - numax, 2) /
                     (2 * np.power(0.25 * numax / cs.FWHM_sigma, 2.0)))
    fitcoef = np.polyfit(xfitdnu, yfitdnu, 1, w=np.sqrt(wfitdnu))
    dnufit = fitcoef[0]
    epsfit = fitcoef[1] / fitcoef[0]
    print(dnufit, epsfit)


def calc_stagger_mlt(temp, grav, metl):
    # Determine mixing-length parameter from 3D simulations
    solarmlt_sbot = 1.98
    solarmlt_dent = 2.09
    x = (temp - 5777.) / 1000.
    y = grav - 4.44

    Fevec = [+0.5, +0.0, -0.5, -1.0, -2.0, -3.0]

    a0vec = (1.97373939, 1.97607815, 1.95635712, 1.96994460, 2.01099658,
             2.13397360)
    a1vec = (-0.13429013, -0.11007132, -0.13364510, -0.14370964, 0.01230849,
             0.05330705)
    a2vec = (0.16320105, 0.17560548, 0.13382454, 0.14900418, 0.16089396,
             0.22228307)
    a3vec = (0.03213152, 0.00397791, 0.02749096, 0.00115447, -0.04127175,
             -0.19292034)
    a4vec = (0.04675864, 0.10333613, 0.04912460, 0.05283657, 0.18048579,
             0.22541197)
    a5vec = (-0.02560510, -0.05869070, -0.04804536, -0.03347076, -0.05957744,
             -0.06493668)
    a6vec = (0.05287061, 0.08055718, 0.05795566, 0.03782256, 0.07440862,
             0.02723039)

    a0int = np.interp(metl, Fevec[::-1], a0vec[::-1])
    a1int = np.interp(metl, Fevec[::-1], a1vec[::-1])
    a2int = np.interp(metl, Fevec[::-1], a2vec[::-1])
    a3int = np.interp(metl, Fevec[::-1], a3vec[::-1])
    a4int = np.interp(metl, Fevec[::-1], a4vec[::-1])
    a5int = np.interp(metl, Fevec[::-1], a5vec[::-1])
    a6int = np.interp(metl, Fevec[::-1], a6vec[::-1])

    mlt_sbot = (a0int +
                (a1int + (a3int + a5int * x + a6int * y) * x + a4int * y) * x +
                a2int * y)

    a0vec = (2.06006479, 2.07706857, 2.08065319, 2.13189578, 2.22904897,
             2.32452703)
    a1vec = (-0.07569695, -0.07928257, -0.11715591, -0.13557833, -0.06863265,
             -0.01166212)
    a2vec = (0.18375048, 0.15337647, 0.13925005, 0.19569437, 0.24814114,
             0.29351541)
    a3vec = (0.01806064, 0.04106182, 0.10587439, 0.03977120, -0.04372929,
             -0.17113568)
    a4vec = (0.16093068, 0.09879519, 0.06301484, 0.10923155, 0.22952285,
             0.30502114)
    a5vec = (-0.11088019, -0.10897204, -0.10459623, -0.07456533, -0.08884621,
             -0.11259495)
    a6vec = (0.16478878, 0.13737732, 0.14323266, 0.11053021, 0.11280482,
             0.07783704)

    a0int = np.interp(metl, Fevec[::-1], a0vec[::-1])
    a1int = np.interp(metl, Fevec[::-1], a1vec[::-1])
    a2int = np.interp(metl, Fevec[::-1], a2vec[::-1])
    a3int = np.interp(metl, Fevec[::-1], a3vec[::-1])
    a4int = np.interp(metl, Fevec[::-1], a4vec[::-1])
    a5int = np.interp(metl, Fevec[::-1], a5vec[::-1])
    a6int = np.interp(metl, Fevec[::-1], a6vec[::-1])

    mlt_dent = (a0int +
                (a1int + (a3int + a5int * x + a6int * y) * x + a4int * y) * x +
                a2int * y)

    print('The computed alpha_mlt from matching the adiabatic entropy value' +
          ' of the deep convection zone is %.2f' % mlt_sbot)
    print('solar alpha_mlt minus this alpha_mlt is %.2f' %
          (solarmlt_sbot - mlt_sbot))
    print('The computed alpha_mlt from matching the entropy jump is %.2f' %
          mlt_dent)
    print('solar alpha_mlt minus this alpha_mlt is %.2f' %
          (solarmlt_dent - mlt_dent))
    return (mlt_sbot, solarmlt_sbot - mlt_sbot,
            mlt_dent, solarmlt_dent - mlt_dent)


def compute_Mv_JC(BT, BT_err, VT, VT_err):
    # https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5
    # Add uncertainties in quadrature as https://arxiv.org/pdf/1510.01731.pdf
    better_unc = np.sqrt(VT_err ** 2 + BT_err ** 2)
    V = [VT - 0.090 * (BT - VT), better_unc]
    print('Johnson magnitude is', V[0], '+/-', V[1])
    return V


def compute_meh(feh, alphaFe):
    falpha = 10 ** alphaFe
    return feh + uncertainties.umath.log10(0.694 * falpha + 0.306)


def compute_asfcorr(evstate, feh, teff, logg, dnu, numax):
    evstate = [evstate]
    feh = [feh]
    teff = [teff]
    logg = [logg]
    dnu = [dnu]
    numax = [numax]

    s = asfgrid.Seism(datadir='/iSIMBA/build/asfgrid_corr/')
    mass, radius = s.get_mass_radius(evstate, feh, teff, dnu, numax)
    logg = s.mr2logg(mass, radius)
    dnu, numax, fdnu = s.get_dnu_numax(evstate, feh, teff, mass,
                                       mass, logg, isfeh=True)
    print('Asf corrected dnu is %.2f µHz ' % dnu[0][0])
    print('Asf corrected numax is %.2f µHz' % numax[0][0])
    print('Asf correction factor is %.2f' % fdnu[0][0])
    return dnu, numax, fdnu


def get_scalingrelation(numax, dnu, teff):
    fnumax = numax / solar.numax
    fdnu = dnu / solar.dnu
    fteff = teff / solar.teff

    r = fnumax * fdnu ** (-2) * fteff ** (1/2)
    m = fnumax ** 3 * fdnu ** (-4) * fteff ** (1.5)
    print('Radius from scaling relation', r)
    print('Mass from scaling relation', m)
    return r, m


def get_mass_from_numaxscalingrelation(r, numax, teff):
    fnumax = numax / solar.numax
    fteff = teff / solar.teff
    m = fnumax * r ** 2 * fteff ** (1/2)
    print('Mass from numax scalingrelation given radius is', m)
    return m


def get_mass_from_dnuscalingrelation(r, dnu):
    fdnu = dnu / solar.dnu
    m = fdnu ** 2 * r ** 3
    print('Mass from dnu scalingrelation given radius is', m)
    return m


def get_sahlholdt2018(numax, dnu, teff):
    fnumax = numax / solar.numax
    fdnu = dnu / solar.dnu
    fteff = teff / solar.teff

    teffcor = teff * 10 ** (-4)
    dnucor = -2.52 * teffcor ** 2 + 2.90 * teffcor + 0.17
    numaxcor = 0.397 * teffcor + 0.771
    r = fnumax * fdnu ** (-2) * fteff ** (1/2)
    m = fnumax ** 3 * fdnu ** (-4) * fteff ** (1.5)
    cr = r * numaxcor ** (-1) * dnucor ** (2)
    cm = m * numaxcor ** (-3) * dnucor ** (4)
    print('Radius from Sahlholdt et al 2018', cr)
    print('Mass from Sahlholdt et al 2018', cm)
    return cr, cm


def get_sharma(numax, dnu, teff, feh, logg, evstate=[0]):
    asfdnu, asfnumax, dnucor = compute_asfcorr(evstate, feh.n, teff.n, logg.n,
                                               dnu.n, numax.n)
    dnucor = dnucor[0][0]
    fnumax = numax / solar.numax
    fdnu = dnu / solar.dnu
    fteff = teff / solar.teff

    r = fnumax * fdnu ** (-2) * dnucor ** (2) * fteff ** (1/2)
    m = fnumax ** 3 * fdnu ** (-4) * dnucor ** (4) * fteff ** (1.5)
    print('Radius from Sharma ', r)
    print('Mass from Sharma', m)
    return r, m


def get_white2011(numax, dnu, teff):
    fnumax = numax / solar.numax
    fdnu = dnu / solar.dnu
    fteff = teff / solar.teff

    kteff = teff * 10 ** (-4)
    teffcorr = -4.29 * kteff ** (2) + 4.84 * kteff - 0.35
    m = teffcorr ** (4) * fnumax ** 3 * fteff ** (3/2) * fdnu ** (-4)
    r = fnumax * fdnu ** (-2) * fteff ** (1/2) * teffcorr ** 2
    print('Mass from White et al 2011', m)
    print('Radius from White et al 2011', r)
    return r, m


def get_surfacegravity(numax, teff):
    fnumax = numax / solar.numax
    fteff = teff / solar.teff

    flogg = (uncertainties.umath.log10(fnumax) +
             0.5 * uncertainties.umath.log10(fteff))
    logg = flogg + solar.logg
    print('The surface gravity is', logg)
    return logg


def get_teff_from_interferometry(fbol, theta):
    stefanboltz = 5.670367 * 1e-5  # erg cm^-2 s^-1 K^4

    teff = ((4 * fbol) / (stefanboltz * theta ** 2)) ** (1/4)
    print('Effective temperature from interferometry is', teff)
    return teff


def get_radius_from_interferometry(plx, theta):
    distance = 1 / plx
    print('The distance is %s parsec' % distance)
    distance *= 4.435 * 1e7  # convert to solar radii
    r = 0.5 * theta * distance
    print('Radius from interferometry is', r)
    return r
