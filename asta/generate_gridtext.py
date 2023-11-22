import os
import h5py
import numpy as np

# Settings (maybe parsing?)
grid = '/home/ADF/stokhoal/BASTA/grids/Garstec_16CygA.hdf5'

assert os.path.exists(grid)

a = h5py.File(grid)

nooftracks = np.ceil(len(a['header/tracks']) / 100) * 100  # round up to nearest 100

# Dimensions
alphaMLT = sorted(np.unique(list(a['header/alphaMLT'])))
alphaMLT = [round(alphaMLT[0], 2), round(alphaMLT[-1], 2)]


# Preamble
preamble = ''
if garstec:
    preamble += r"\newcommand{\Garstec}{\textsc{Garstec}\xspace}"
if adipls:
    preamble += r"\newcommand{\adipls}{\textsc{ADIPLS}\xspace}"
preamble += r"\newcommand{\range}[2]{\ensuremath{\ni ] #1;#2 [}\xspace}"

# Paragraph
paragraph = f" We determine the stellar properties using grid-based stellar modelling."
paragraph += "We construct grids of theoretical models of stellar evolution covering the necessary parameter space, which we compare the observed parameters to the predicted theoretical quantities. Quantities like stellar age can thus be inferred from the constraints given by the other observables."


# Sampling
text = ''
text += f'We build a grid with ${\sim}{nooftracks}$ evolutionary tracks of stellar models'
if sobol:
    text += f"We sampled the parameter space by utilising the Sobol quasi-random, low-discrepancy sequences to uniformly populate the parameter space \citep{sobol1,sobol2,sobol6,sobol4,sobol5,sobol3}."

# Stellar evolution code
if garstec:
    text += f"Grids of stellar models were computed with the Garching Stellar Evolution Code \citep[\Garstec][]{weiss2008}."
else:
    raise ValueError, 'What is your stellar evolution code?'


# EOS
text += r"The code utilises a combination of the equation of state by the \textsc{OPAL} group \citep{rogers1996,rogers2002}, and the Mihalas-Hummer-Däppen equation of state \citep{mihalas1988,hummer1988,daeppen1988,mihalas1990}."

# Opacities
text += "We use the OPAL opacities \citep{rogers1992,iglesias1996} at high temperatures supplemented by the opacities of \citet{ferguson2005} at low temperatures."
text += r"\Garstec uses the NACRE nuclear reaction rates \citep{angulo1999} except for $^{14}$N($p,\gamma$)$^{15}$O and $^{12}$C($\alpha,\gamma$)$^{16}$O for which the rates from \citet{formicola2004} and \citet{hammer2005} were used."

# Solar mixture and opacities
if GV98:
    text += f'The stellar models are computed using the \citet{grevesse1998} solar mixture'
elif ASP09:
    text += f'The stellar models are computed using the \citet{asplund2009} solar mixture'
else:
    raise ValueError, 'What is your solar mixture?'

# Convection
text += f"Convection in the models are parameterised using mixing-length theory \citep{bohm1958, kippenhahn2012},"
if convection:
    text += f"where the mixing length parameter is allowed to vary in the range \range[{alphaMLT[0]}][{alphaMLT[1]}]."
else: 
    if len(np.unique(alphaMLT)) == 1:
       text += f"where the mixing length parameter is kept constant at the solar-calibrated value of $\amlt={alphaMLT[0]}$ as determined by a standard solar model calibration."

# ODEA
if normal:
    text += "Diffusion and settling of helium and heavier elements were not included, neither was convective overshooting."

assert normal or (dif or ove or massloss or alphaFe)
if massloss:
    text += f"The mass loss ($\eta$) ranges from $0.0$ to $0.3$ following the \citet{reimers1977} formalism."
if dif:
    text += r"Atomic diffusion of elements are treated following the prescription by \cite{thoul1994}."
if ove:

if alphaFe:
    text += r''

if chemevo:
    text = r'The primordial helium is assumed to be 0.248 \citep{fields2020} and the Galactic chemical enrichment law is assumed to be $\Delta Y/\Delta Z=1.4$.'


# Atmosphere
text += "We used an Eddington grey atmosphere."

# Dimensions
The stellar grid samples masses from $1.0$--$1.5$~\si{\solarmass} in steps of \SI{0.01}{\solarmass} and it samples metallicities from $\meh=-0.32$ to $-0.14$ in steps of $0.03$, assuming a fixed linear Galactic chemical evolution model of $\Delta Y / \Delta Z = 1.4$ \citep{balser2006}.

# Dnu range
The grid covers \dnu in the range $50$--$60$~\si{\micro\hertz}, thus spanning the parameter space from about $\SI{5}{\micro\hertz}$ on both sides of the observed \dnu.

# Individual frequencies
if adipls:
    text += r'The theoretical oscillation frequencies of the models are computed using the Aarhus adiabatic oscillation package \adipls \citep[][]{jcd2008}.'

if parallax:
    text += r'For the computation of synthetic magnitudes, we use the bolometric corrections of \citet{hidalgo2018}'

# BASTA
text += r'We use the BAyesian STellar Algorithm \citep[\basta;][]{aguirresilvaaguirre2015,silvaaguirre2017} to determine the stellar parameters. Given a precomputed grid of stellar models, \basta uses a Bayesian approach to compute the posterior distribution of a given stellar parameter using a set of observational constraints. '

if imfprior:
    text = r'\basta allows for prior probability distributions to be taken into account when computing the posterior distributions. We used the Salpeter initial mass function \citep{salpeter1955} to quantify the expected mass distribution of stars favouring low-mass stars as the most abundant.'

if universeageprior:  # This is directly from Borre et al. 2021
    text += r'Additionally, we include an upper limit on the stellar ages of 15 Gyr. This is done to avoid nonphysical solutions for stars older than the age of the Universe. Despite the solutions not being physical at above the age of the Universe (13.7 Gyr), they can still hold statistical significance and we do, therefore, not truncate the solutions at 13.7 Gyr but allow them to stretch to 15 Gyr. For the remaining parameters we use uniform priors.'


# Fitting


# Bibliography
bib = ''
if sobol:
    bib += """
    @ARTICLE{sobol1,
        author = {{Sobol}, Ilya and {Levithan}, Y.L.},
        title = "{The Production of Points Uniformly Distributed in a Multidimensional Cube (in Russian)}",
        journal = {IPM Akademii Nauk SSSR},
        number = {40},
        year = {1976}
    }

    @ARTICLE{sobol2,
        author = {{Sobol}, Ilya},
        title = "{Uniformly Distributed Sequences with an Additional Uniform Property}",
        journal = {USSR Computational Mathematics and Mathematical Physics},
        volume = {16},
        year = {1977},
        pages = {236-242}
    }

    @ARTICLE{sobol3,
        author = {{Joe}, Stephen and {Kuo}, Frances},
        title = "{Remark on Algorithm 659: Implementing Sobol's Quasirandom Sequence Generator}",
        journal = {ACM Transactions on Mathematical Software},
        number = {1},
        volume = {29},
        year = {2003},
        month = {3},
        pages = {49-57}
    }

    @ARTICLE{sobol4,
        author = {{Fox}, Bennett},
        title = "{Algorithm 647: Implementation and Relative Efficiency of Quasirandom Sequence Generators}",
        journal = {ACM Transactions on Mathematical Software},
        number = {4},
        volume = {12},
        year = {1986},
        month = {12},
        pages = {362-376}
    }

    @ARTICLE{sobol5,
        author = {{Bratley}, Paul and {Fox}, Bennett},
        title = "{Algorithm 659: Implementing Sobol's Quasirandom Sequence Generator}",
        journal = {ACM Transactions on Mathematical Software},
        number = {1},
        volume = {14},
        year = {1988},
        month = {3},
        pages = {88-100}
    }

    @ARTICLE{sobol6,
        author = {{Antonov}, I.A. and {Saleev}, V.M.},
        title = "{An Economic Method of Computing LP Tau-Sequences}",
        journal = {USSR Computational Mathematics and Mathematical Physics},
        volume = {19},
        year = {1980},
        pages = {252-256}
    }
"""

if GV98:
    bib += """
    """
elif ASP09:
    bib += r"""
    @ARTICLE{asplund2009,
       author = {{Asplund}, Martin and {Grevesse}, Nicolas and {Sauval}, A. Jacques and {Scott}, Pat},
        title = "{The Chemical Composition of the Sun}",
      journal = {\araa},
     keywords = {Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2009,
        month = sep,
       volume = {47},
       number = {1},
        pages = {481-522},
          doi = {10.1146/annurev.astro.46.060407.145222},
archivePrefix = {arXiv},
       eprint = {0909.0948},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2009ARA&A..47..481A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    """

if garstec:
    bib += r"""
    @article{weiss2008,
      author = {{Weiss}, A. and {Schlattl}, H.},
    title = "{GARSTEC: the Garching Stellar Evolution Code. The direct descendant of the legendary Kippenhahn code}",
    journal = {￼￼apss},
    year = 2008,
    month = 8,
    volume = 316,
    pages = {99-106},
    doi = {10.1007/s10509-007-9606-5},
    adsurl = {http://adsabs.harvard.edu/abs/2008Ap26SS.316...99W},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    """

bib += r"""
@ARTICLE{rogers1996,
	author = {{Rogers}, Forrest J. and {Swenson}, Fritz J. and {Iglesias}, Carlos A.},
	title = "{OPAL Equation-of-State Tables for Astrophysical Applications}",
	journal = {\apj},
	keywords = {ATOMIC PROCESSES, EQUATION OF STATE, ATOMIC DATA},
	year = 1996,
	month = jan,
	volume = {456},
	pages = {902},
	doi = {10.1086/176705},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1996ApJ...456..902R},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{rogers2002,
	author = {{Rogers}, F.~J. and {Nayfonov}, A.},
	title = "{Updated and Expanded OPAL Equation-of-State Tables: Implications for Helioseismology}",
	journal = {\apj},
	keywords = {Atomic Processes, Equation of State, Sun: Oscillations},
	year = 2002,
	month = sep,
	volume = {576},
	number = {2},
	pages = {1064-1074},
	doi = {10.1086/341894},
	adsurl = {https://ui.adsabs.harvard.edu/abs/2002ApJ...576.1064R},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{daeppen1988,
	author = {{Daeppen}, Werner and {Mihalas}, Dimitri and {Hummer}, D.~G. and
	{Mihalas}, Barbara Weibel},
	title = "{The Equation of State for Stellar Envelopes. III. Thermodynamic Quantities}",
	journal = {\apj},
	keywords = {Equations Of State, Free Energy, Helium Plasma, Hydrogen Plasma, Stellar Envelopes, Thermodynamic Properties, Helium Hydrogen Atmospheres, Internal Energy, Ionized Gases, Isotherms, Partial Differential Equations, Astrophysics, ATOMIC PROCESSES, EQUATION OF STATE, STARS: INTERIORS},
	year = 1988,
	month = sep,
	volume = {332},
	pages = {261},
	doi = {10.1086/166650},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1988ApJ...332..261D},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{hummer1988,
	author = {{Hummer}, D.~G. and {Mihalas}, Dimitri},
	title = "{The Equation of State for Stellar Envelopes. I. an Occupation Probability Formalism for the Truncation of Internal Partition Functions}",
	journal = {\apj},
	keywords = {Equations Of State, Stellar Atmospheres, Stellar Envelopes, Stellar Interiors, Charged Particles, Cosmic Plasma, Coulomb Potential, Perturbation Theory, Probability Density Functions, Astrophysics, ATOMIC PROCESSES, EQUATION OF STATE, STARS: ATMOSPHERES},
	year = 1988,
	month = aug,
	volume = {331},
	pages = {794},
	doi = {10.1086/166600},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1988ApJ...331..794H},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{mihalas1988,
	author = {{Mihalas}, Dimitri and {Dappen}, Werner and {Hummer}, D.~G.},
	title = "{The Equation of State for Stellar Envelopes. II. Algorithm and Selected Results}",
	journal = {\apj},
	keywords = {Charged Particles, Computational Astrophysics, Equations Of State, Stellar Atmospheres, Stellar Envelopes, Algorithms, Cosmic Plasma, Coulomb Collisions, Free Energy, Perturbation Theory, Astrophysics, ATOMIC PROCESSES, EQUATION OF STATE, STARS: ATMOSPHERES},
	year = 1988,
	month = aug,
	volume = {331},
	pages = {815},
	doi = {10.1086/166601},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1988ApJ...331..815M},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{mihalas1990,
	author = {{Mihalas}, Dimitri and {Hummer}, D.~G. and {Mihalas}, Barbara Weibel and
	{Daeppen}, Werner},
	title = "{The Equation of State for Stellar Envelopes. IV. Thermodynamic Quantities and Selected Ionization Fractions for Six Elemental Mixes}",
	journal = {\apj},
	keywords = {Equations Of State, Gas Ionization, Metallic Stars, Stellar Envelopes, Thermodynamics, Helium Ions, Hydrogen Ions, Temperature Dependence, Astrophysics, ATOMIC PROCESSES, EQUATION OF STATE, STARS: ATMOSPHERES},
	year = 1990,
	month = feb,
	volume = {350},
	pages = {300},
	doi = {10.1086/168383},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1990ApJ...350..300M},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{rogers1992,
	author = {{Rogers}, Forrest J. and {Iglesias}, Carlos A.},
	title = "{Radiative Atomic Rosseland Mean Opacity Tables}",
	journal = {\apjs},
	keywords = {Atomic Spectra, Opacity, Radiative Transfer, Stellar Composition, Stellar Envelopes, Stellar Interiors, Abundance, Equations Of State, Metallicity, Stellar Atmospheres, Atomic and Molecular Physics, ATOMIC DATA, ATOMIC PROCESSES, STARS: INTERIORS},
	year = 1992,
	month = apr,
	volume = {79},
	pages = {507},
	doi = {10.1086/191659},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1992ApJS...79..507R},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{iglesias1996,
	author = {{Iglesias}, Carlos A. and {Rogers}, Forrest J.},
	title = "{Updated Opal Opacities}",
	journal = {\apj},
	keywords = {ATOMIC DATA, ATOMIC PROCESSES, STARS: INTERIORS},
	year = 1996,
	month = jun,
	volume = {464},
	pages = {943},
	doi = {10.1086/177381},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1996ApJ...464..943I},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{ferguson2005,
	author = {{Ferguson}, Jason W. and {Alexander}, David R. and {Allard}, France and
	{Barman}, Travis and {Bodnarik}, Julia G. and {Hauschildt}, Peter H. and
	{Heffner-Wong}, Amanda and {Tamanai}, Akemi},
	title = "{Low-Temperature Opacities}",
	journal = {\apj},
	keywords = {Atomic Data, Equation of State, Methods: Numerical, Molecular Data, Astrophysics},
	year = 2005,
	month = apr,
	volume = {623},
	number = {1},
	pages = {585-596},
	doi = {10.1086/428642},
	archivePrefix = {arXiv},
	eprint = {astro-ph/0502045},
	primaryClass = {astro-ph},
	adsurl = {https://ui.adsabs.harvard.edu/abs/2005ApJ...623..585F},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{angulo1999,
	author = {{Angulo}, C. and {Arnould}, M. and {Rayet}, M. and {Descouvemont}, P. and
	{Baye}, D. and {Leclercq-Willain}, C. and {Coc}, A. and {Barhoumi}, S. and
	{Aguer}, P. and {Rolfs}, C. and {Kunz}, R. and {Hammer}, J.~W. and
	{Mayer}, A. and {Paradellis}, T. and {Kossionides}, S. and
	{Chronidou}, C. and {Spyrou}, K. and {degl'Innocenti}, S. and
	{Fiorentini}, G. and {Ricci}, B. and {Zavatarelli}, S. and
	{Providencia}, C. and {Wolters}, H. and {Soares}, J. and {Grama}, C. and
	{Rahighi}, J. and {Shotter}, A. and {Lamehi Rachti}, M.},
	title = "{A compilation of charged-particle induced thermonuclear reaction rates}",
	journal = {\nphysa},
	year = 1999,
	month = aug,
	volume = {656},
	number = {1},
	pages = {3-183},
	doi = {10.1016/S0375-9474(99)00030-5},
	adsurl = {https://ui.adsabs.harvard.edu/abs/1999NuPhA.656....3A},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@article{formicola2004,
    title = {Astrophysical S-factor of 14N(p,γ)15O},
    journal = {Physics Letters B},
    volume = {591},
    number = {1},
    pages = {61-68},
    year = {2004},
    issn = {0370-2693},
    doi = {https://doi.org/10.1016/j.physletb.2004.03.092},
    url = {https://www.sciencedirect.com/science/article/pii/S037026930400601X},
    author = {A Formicola and G Imbriani and H Costantini and C Angulo and D Bemmerer and R Bonetti and C Broggini and P Corvisiero and J Cruz and P Descouvemont and Z Fülöp and G Gervino and A Guglielmetti and C Gustavino and G Gyürky and A.P Jesus and M Junker and A Lemut and R Menegazzo and P Prati and V Roca and C Rolfs and M Romano and C {Rossi Alvarez} and F Schümann and E Somorjai and O Straniero and F Strieder and F Terrasi and H.P Trautvetter and A Vomiero and S Zavatarelli}
    }

@article{hammer2005,
    title = {E1 and E2 capture cross section and astrophysical reaction rate of the key reaction 12C(α,γ)16O},
    journal = {Nuclear Physics A},
    volume = {758},
    pages = {363-366},
    year = {2005},
    note = {Nuclei in the Cosmos VIII},
    issn = {0375-9474},
    doi = {https://doi.org/10.1016/j.nuclphysa.2005.05.066},
    url = {https://www.sciencedirect.com/science/article/pii/S0375947405007153},
    author = {J.W. Hammer and M. Fey and R. Kunz and J. Kiener and V. Tatischeff and F. Haas and J.L. Weil and M. Assunção and C. Beck and C. Boukari-Pelissie and A. Coc and J.J. Correia and S. Courtin and F. Fleurot and E. Galanopoulos and C. Grama and F. Hammache and S. Harissopulos and A. Korichi and E. Krmpotić and D. {Le Du} and A. Lopez-Martens and D. Malcherek and R. Meunier and P. Papka and T. Paradellis and M. Rousseau and N. Rowley and G. Staudt and S. Szilner}
    }

￼
@ARTICLE{bohm1958,
       author = {{B{\"o}hm-Vitense}, E.},
        title = "{{\"U}ber die Wasserstoffkonvektionszone in Sternen verschiedener Effektivtemperaturen und Leuchtkr{\"a}fte. Mit 5 Textabbildungen}",
      journal = {\zap},
         year = 1958,
        month = jan,
       volume = {46},
        pages = {108},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1958ZA.....46..108B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

￼
@BOOK{kippenhahn2012,
       author = {{Kippenhahn}, Rudolf and {Weigert}, Alfred and {Weiss}, Achim},
        title = "{Stellar Structure and Evolution}",
         year = 2012,
          doi = {10.1007/978-3-642-30304-3},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2012sse..book.....K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{reimers1977,
       author = {{Reimers}, D.},
        title = "{Observational evidence for mass-loss from K giants, G and K supergiants.}",
      journal = {\aap},
     keywords = {G Stars, Giant Stars, K Stars, Red Giant Stars, Stellar Mass Ejection, Stellar Winds, Supergiant Stars, Absorption Spectra, Line Spectra, M Stars, Stellar Envelopes, Stellar Spectra, Astrophysics},
         year = 1977,
        month = may,
       volume = {57},
        pages = {395-400},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1977A&A....57..395R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{thoul1994,
       author = {{Thoul}, Anne A. and {Bahcall}, John N. and {Loeb}, Abraham},
        title = "{Element Diffusion in the Solar Interior}",
      journal = {\apj},
     keywords = {Abundance, Diffusion, Flow Equations, Heavy Elements, Helium, Solar Interior, Stellar Composition, Stellar Models, Computerized Simulation, Stellar Evolution, Subroutines, Solar Physics, DIFFUSION, STARS: ABUNDANCES, STARS: INTERIORS, SUN: INTERIOR, Astrophysics},
         year = 1994,
        month = feb,
       volume = {421},
        pages = {828},
          doi = {10.1086/173695},
archivePrefix = {arXiv},
       eprint = {astro-ph/9304005},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1994ApJ...421..828T},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

￼
@ARTICLE{jcd2008,
	author = {{Christensen-Dalsgaard}, J{\o}rgen},
	title = "{ADIPLS{\textemdash}the Aarhus adiabatic oscillation package}",
	journal = {\apss},
	keywords = {Stars: oscillations, Numerical methods, Asteroseismology, Astrophysics},
	year = 2008,
	month = aug,
	volume = {316},
	number = {1-4},
	pages = {113-120},
	doi = {10.1007/s10509-007-9689-z},
	archivePrefix = {arXiv},
	eprint = {0710.3106},
	primaryClass = {astro-ph},
	adsurl = {https://ui.adsabs.harvard.edu/abs/2008Ap&SS.316..113C},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
￼
@ARTICLE{ball2014,
	author = {{Ball}, W.~H. and {Gizon}, L.},
	title = "{A new correction of stellar oscillation frequencies for near-surface effects}",
	journal = {\aap},
	keywords = {asteroseismology, stars: oscillations, stars: individual: HD 52265, Astrophysics - Solar and Stellar Astrophysics},
	year = 2014,
	month = aug,
	volume = {568},
	eid = {A123},
	pages = {A123},
	doi = {10.1051/0004-6361/201424325},
	archivePrefix = {arXiv},
	eprint = {1408.0986},
	primaryClass = {astro-ph.SR},
	adsurl = {https://ui.adsabs.harvard.edu/abs/2014A&A...568A.123B},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

"""
print('WARNING: DO NOT USE THIS TEXT FOR PAPERS ETC AS IT IS BASICALLY COPY-PASTE')
print('PLEASE DOUBLE-CHECK IF THE TEXT MATCHES YOUR EXPECTATIONS')
print('TO BE USED FOR EASILY COMMUNICATING A SET-UP PRE-PUBLICATION')
print(preamble)
print('')
print(paragraph)
print(text)
print('')
print('')
print(bib)
