import os
import h5py
import numpy as np

# Settings (maybe parsing?)
grid = "/home/ADF/stokhoal/BASTA/grids/Garstec_16CygA.hdf5"
outdir = "./gridtext/"
if not os.path.exists(outdir):
    os.path.makedirs(outdir)
filename = os.path.join(outdir, f"griddescription_{grid.split('/')[-1][:-5]}.txt")
bibname = os.path.join(outdir, f"bib_{grid.split('/')[-1][:-5]}.bib")
sobol = True
mindnu = None
maxdnu = None
parallax = True
imfprior = True
universeageprior = False

assert os.path.exists(grid)

a = h5py.File(grid)
nooftracks = np.ceil(len(a["header/tracks"]) / 100) * 100  # round up to nearest 100


# Preamble
preamble = r"""
\usepackage{siunitx}    % Typesetting of numbers and/or units
\usepackage{xspace}     % Appropriate space for commans
\DeclareSIUnit{\solarmass}{\ensuremath{M_\odot}\xspace}

\newcommand{\basta}{\textsc{basta}\xspace}
"""
if "Garstec" in grid:
    preamble += r"\newcommand{\Garstec}{\textsc{Garstec}\xspace}"
    preamble += r"\newcommand{\adipls}{\textsc{adipls}\xspace}"
preamble += r"""
\newcommand{\dnu}{\ensuremath{\Delta\nu}\xspace}
\newcommand{\numax}{\ensuremath{\nu_\textup{max}}\xspace}
\newcommand{\feh}{\ensuremath{[\textup{Fe}/\textup{H}]}\xspace}
\newcommand{\meh}{\ensuremath{[\textup{M}/\textup{H}]}\xspace}
\newcommand{\teff}{\ensuremath{\var{T}{eff}}\xspace}
\newcommand{\logg}{\ensuremath{\log g}\xspace}
\newcommand{\lphot}{\ensuremath{\var{L}{Phot}}\xspace}
"""

# Paragraph
paragraph = r"We determine the stellar properties using grid-based stellar modelling. "
paragraph += "We construct grids of theoretical models of stellar evolution covering the necessary parameter space, which we compare the observed parameters to the predicted theoretical quantities. Quantities like stellar age can thus be inferred from the constraints given by the other observables."


# Sampling
text = (
    f"We build a grid with $\sim{nooftracks}$ evolutionary tracks of stellar models. "
)
if sobol:
    text += r"We sampled the parameter space by utilising the Sobol quasi-random, low-discrepancy sequences to uniformly populate the parameter space \citep{sobol1,sobol2,sobol6,sobol4,sobol5,sobol3}. "

# Stellar evolution code
if "Garstec" in grid:
    text += r"Grids of stellar models were computed with the Garching Stellar Evolution Code \citep[\Garstec][]{weiss2008}. "
else:
    raise ValueError("What is your stellar evolution code?")


# EOS
text += r"The code utilises a combination of the equation of state by the \textsc{OPAL} group \citep{rogers1996,rogers2002}, and the Mihalas-Hummer-Däppen equation of state \citep{mihalas1988,hummer1988,daeppen1988,mihalas1990}. "

# Opacities
text += "We use the OPAL opacities \citep{rogers1992,iglesias1996} at high temperatures supplemented by the opacities of \citet{ferguson2005} at low temperatures. "
text += r"\Garstec uses the NACRE nuclear reaction rates \citep{angulo1999} except for $^{14}$N($p,\gamma$)$^{15}$O and $^{12}$C($\alpha,\gamma$)$^{16}$O for which the rates from \citet{formicola2004} and \citet{hammer2005} were used. "

# Solar mixture and opacities
if "_gv98_" in list(a["solar_models"])[0]:
    text += r"The stellar models are computed using the \citet{grevesse1998} solar mixture. "
elif "_as09_" in list(a["solar_models"])[0]:
    text += (
        r"The stellar models are computed using the \citet{asplund2009} solar mixture. "
    )
else:
    raise ValueError("What is your solar mixture?")

text += "\n"

# Convection
text += r"Convection in the models are parameterised using mixing-length theory \citep{bohm1958, kippenhahn2012},"
alphaMLT = sorted(np.unique(list(a["header/alphaMLT"])))
alphaMLT = [round(alphaMLT[0], 2), round(alphaMLT[-1], 2)]
if len(np.unique(alphaMLT)) > 1:
    text += f"where the mixing length parameter is allowed to vary in the range {alphaMLT[0]}--{alphaMLT[1]}. "
else:
    text += r"where the mixing length parameter is kept constant at the solar-calibrated value of $\amlt="
    text += f"{alphaMLT[0]}$ as determined by a standard solar model calibration. "

# ODEA
if {"dif", "ove"} <= set(a["header/pars_constant"]):
    text += "Diffusion and settling of helium and heavier elements were not included, neither was convective overshooting. "

text = r"The primordial helium is assumed to be $0.248$ \citep{fields2020} and we use a fixed the Galactic chemical enrichment law of $\Delta Y/\Delta Z=1.4$ \citep{balser2006}. "

# Atmosphere
text += "We used an Eddington grey atmosphere. "

# Dimensions
params = {
    "massini": ["masses", "\si{\solarmass}"],
    "FeHini": ["initial iron abundances", "dex"],
    "MeHini": ["initial metallicity", "dex"],
    "alphaMLT": ["mixing-length parameter", ""],
    "yini": ["initial helium fraction", ""],
    "alphaFe": [r"$\alpha$ enhancement", ""],
    "dif": [],
    "eta": ["mass loss efficiency", " following the \citet{reimers1977} formalism."],
    "ove": ["overshooting efficiency", ""],
}
text += f"The stellar grid samples "
for i, param in enumerate(list(a["header/pars_sampled"])):
    param = param.astype(str)
    if param == "dif":
        text += r"Atomic diffusion of elements are treated following the prescription by \cite{thoul1994}."
        continue
    paramlist = sorted(np.unique(list(a[f"header/{param}"])))
    paramlist = [round(paramlist[0], 3), round(paramlist[-1], 3)]
    text += (
        f"{params[param][0]} from ${paramlist[0]}$--${paramlist[1]}$~{params[param][1]}"
    )
    if i == len(list(a["header/pars_sampled"])):
        text += ". "
    else:
        text += ", "

# Dnu range
if mindnu is not None:
    text += f"The grid covers \dnu in the range ${mindnu}$--${maxdnu}$"
    text += r"~\si{\micro\hertz}. "
elif maxdnu is not None:
    text += r"The range of models in a given track is limited between the zero-age main-sequence and when the model reaches a large frequency separation of $\dnu = {maxdnu}~\si{\micro\hertz}$. Here, the zero-age main-sequence is determined as where the ratio between the hydrogen burning luminosity and the total luminosity reaches $1$. "

# Individual frequencies
if "Garstec" in grid:
    text += r"The theoretical oscillation frequencies of the models are computed using the Aarhus adiabatic oscillation package \citep[\adipls;][]{jcd2008}. "

if parallax:
    text += r"For the computation of synthetic magnitudes, we use the bolometric corrections of \citet{hidalgo2018} and we use the dust map from \citet{green2019} for computing the extinction. "

text += "\n"

# BASTA
text += r"We use the BAyesian STellar Algorithm \citep[\basta;][]{silvaaguirre2015,silvaaguirre2017,aguirreborsenkoch2022} to determine the stellar parameters. Given a precomputed grid of stellar models, \basta uses a Bayesian approach to compute the posterior distribution of a given stellar parameter using a set of observational constraints. "

if imfprior:
    text += r"\basta allows for prior probability distributions to be taken into account when computing the posterior distributions. We used the Salpeter initial mass function \citep{salpeter1955} to quantify the expected mass distribution of stars favouring low-mass stars as the most abundant. "

if universeageprior:  # This is directly from Borre et al. 2021
    text += r"Additionally, we include an upper limit on the stellar ages of 15 Gyr. This is done to avoid nonphysical solutions for stars older than the age of the Universe. Despite the solutions not being physical at above the age of the Universe (13.7 Gyr), they can still hold statistical significance and we do, therefore, not truncate the solutions at 13.7 Gyr but allow them to stretch to 15 Gyr. For the remaining parameters we use uniform priors. "


# Bibliography
bib = ""
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

if "_gv98_" in list(a["solar_models"])[0]:
    bib += r"""
@article{grevesse1998,
    author = {{Grevesse}, N. and {Sauval}, A.~J.},
    title = "{Standard Solar Composition}",
    journal = {\ssr},
    keywords = {Sun: abundances, Meteorites: abundances, Solar spectroscopy},
    year = 1998,
    month = may,
    volume = 85,
    pages = {161-174},
    doi = {10.1023/A:1005161325181},
    adsurl = {http://adsabs.harvard.edu/abs/1998SSRv...85..161G},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
        """
elif "_as09_" in list(a["solar_models"])[0]:
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

if "Garstec" in grid:
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


@article{hidalgo2018,
	doi = {10.3847/1538-4357/aab158},
	url = {https://doi.org/10.3847/1538-4357/aab158},
	year = 2018,
	month = {mar},
	publisher = {American Astronomical Society},
	volume = {856},
	number = {2},
	pages = {125},
	author = {Sebastian L. Hidalgo and Adriano Pietrinferni and Santi Cassisi and Maurizio Salaris and Alessio Mucciarelli and Alessandro Savino and Antonio Aparicio and Victor Silva Aguirre and Kuldeep Verma},
	title = {The Updated {BaSTI} Stellar Evolution Models and Isochrones. I. Solar-scaled Calculations},
	journal = {The Astrophysical Journal},
}


@ARTICLE{green2015,
       author = {{Green}, Gregory M. and {Schlafly}, Edward F. and
         {Finkbeiner}, Douglas P. and {Rix}, Hans-Walter and {Martin}, Nicolas and
         {Burgett}, William and {Draper}, Peter W. and {Flewelling}, Heather and
         {Hodapp}, Klaus and {Kaiser}, Nicholas and {Kudritzki}, Rolf Peter and
         {Magnier}, Eugene and {Metcalfe}, Nigel and {Price}, Paul and
         {Tonry}, John and {Wainscoat}, Richard},
        title = "{A Three-dimensional Map of Milky Way Dust}",
      journal = {\apj},
     keywords = {dust, extinction, Galaxy: structure, methods: statistical, Astrophysics - Astrophysics of Galaxies},
         year = 2015,
        month = sep,
       volume = {810},
       number = {1},
          eid = {25},
        pages = {25},
          doi = {10.1088/0004-637X/810/1/25},
archivePrefix = {arXiv},
       eprint = {1507.01005},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015ApJ...810...25G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{green2019,
       author = {{Green}, Gregory M. and {Schlafly}, Edward and {Zucker}, Catherine and
         {Speagle}, Joshua S. and {Finkbeiner}, Douglas},
        title = "{A 3D Dust Map Based on Gaia, Pan-STARRS 1, and 2MASS}",
      journal = {\apj},
     keywords = {Interstellar reddening, Interstellar dust extinction, Galaxy structure, Galaxy stellar content, Interstellar dust, 853, 837, 622, 621, 836, Astrophysics - Astrophysics of Galaxies},
         year = 2019,
        month = dec,
       volume = {887},
       number = {1},
          eid = {93},
        pages = {93},
          doi = {10.3847/1538-4357/ab5362},
archivePrefix = {arXiv},
       eprint = {1905.02734},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...887...93G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{aguirre2015,
       author = {{Silva Aguirre}, V. and {Davies}, G.~R. and {Basu}, S. and
         {Christensen-Dalsgaard}, J. and {Creevey}, O. and {Metcalfe}, T.~S. and
         {Bedding}, T.~R. and {Casagrande}, L. and {Handberg}, R. and
         {Lund}, M.~N. and {Nissen}, P.~E. and {Chaplin}, W.~J. and {Huber}, D. and
         {Serenelli}, A.~M. and {Stello}, D. and {Van Eylen}, V. and
         {Campante}, T.~L. and {Elsworth}, Y. and {Gilliland}, R.~L. and
         {Hekker}, S. and {Karoff}, C. and {Kawaler}, S.~D. and {Kjeldsen}, H. and
         {Lundkvist}, M.~S.},
        title = "{Ages and fundamental properties of Kepler exoplanet host stars from asteroseismology}",
      journal = {\mnras},
     keywords = {asteroseismology, planets and satellites: fundamental parameters, stars: evolution, stars: fundamental parameters, stars: oscillations, planetary systems, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2015,
        month = sep,
       volume = {452},
       number = {2},
        pages = {2127-2148},
          doi = {10.1093/mnras/stv1388},
archivePrefix = {arXiv},
       eprint = {1504.07992},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2127S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{aguirre2017,
       author = {{Silva Aguirre}, V. and {Lund}, Mikkel N. and {Antia}, H.~M. and
         {Ball}, Warrick H. and {Basu}, Sarbani and
         {Christensen-Dalsgaard}, J{\o}rgen and {Lebreton}, Yveline and
         {Reese}, Daniel R. and {Verma}, Kuldeep and {Casagrande}, Luca and
         {Justesen}, Anders B. and {Mosumgaard}, Jakob R. and
         {Chaplin}, William J. and {Bedding}, Timothy R. and {Davies}, Guy R. and
         {Handberg}, Rasmus and {Houdek}, G{\"u}nter and {Huber}, Daniel and
         {Kjeldsen}, Hans and {Latham}, David W. and {White}, Timothy R. and
         {Coelho}, Hugo R. and {Miglio}, Andrea and {Rendle}, Ben},
        title = "{Standing on the Shoulders of Dwarfs: the Kepler Asteroseismic LEGACY Sample. II.Radii, Masses, and Ages}",
      journal = {\apj},
     keywords = {asteroseismology, stars: fundamental parameters, stars: oscillations, Astrophysics - Solar and Stellar Astrophysics},
         year = 2017,
        month = feb,
       volume = {835},
       number = {2},
          eid = {173},
        pages = {173},
          doi = {10.3847/1538-4357/835/2/173},
archivePrefix = {arXiv},
       eprint = {1611.08776},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2017ApJ...835..173S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{aguirreborsenkoch2022,
       author = {{Aguirre B{\o}rsen-Koch}, V. and {R{\o}rsted}, J.~L. and {Justesen}, A.~B. and {Stokholm}, A. and {Verma}, K. and {Winther}, M.~L. and {Knudstrup}, E. and {Nielsen}, K.~B. and {Sahlholdt}, C. and {Larsen}, J.~R. and {Cassisi}, S. and {Serenelli}, A.~M. and {Casagrande}, L. and {Christensen-Dalsgaard}, J. and {Davies}, G.~R. and {Ferguson}, J.~W. and {Lund}, M.~N. and {Weiss}, A. and {White}, T.~R.},
        title = "{The BAyesian STellar algorithm (BASTA): a fitting tool for stellar studies, asteroseismology, exoplanets, and Galactic archaeology}",
      journal = {\mnras},
     keywords = {asteroseismology, methods: numerical, methods: statistical, stars: fundamental parameters, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2022,
        month = jan,
       volume = {509},
       number = {3},
        pages = {4344-4364},
          doi = {10.1093/mnras/stab2911},
archivePrefix = {arXiv},
       eprint = {2109.14622},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.4344A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


"""
if imfprior:
    bib += r"""
@ARTICLE{salpeter1955,
       author = {{Salpeter}, Edwin E.},
        title = "{The Luminosity Function and Stellar Evolution.}",
      journal = {\apj},
         year = 1955,
        month = jan,
       volume = {121},
        pages = {161},
          doi = {10.1086/145971},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
bib += r"""
@ARTICLE{balser2006,
   author = {{Balser}, D.~S.},
    title = "{The Chemical Evolution of Helium}",
  journal = {\aj},
   eprint = {astro-ph/0608436},
 keywords = {ISM: H II Regions, ISM: Abundances, Radio Lines: ISM},
     year = 2006,
    month = dec,
   volume = 132,
    pages = {2326-2332},
      doi = {10.1086/508515},
   adsurl = {http://adsabs.harvard.edu/abs/2006AJ....132.2326B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
bib += r"""
@ARTICLE{fields2020,
       author = {{Fields}, Brian D. and {Olive}, Keith A. and {Yeh}, Tsung-Han and {Young}, Charles},
        title = "{Big-Bang Nucleosynthesis after Planck}",
      journal = {\jcap},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, High Energy Physics - Phenomenology, Nuclear Experiment},
         year = 2020,
        month = mar,
       volume = {2020},
       number = {3},
          eid = {010},
        pages = {010},
          doi = {10.1088/1475-7516/2020/03/010},
archivePrefix = {arXiv},
       eprint = {1912.01132},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020JCAP...03..010F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

with open(filename, "w") as outfile:
    outfile.write(
        "WARNING: DO NOT USE THIS TEXT FOR PAPERS ETC AS IT IS BASICALLY COPY-PASTE"
    )
    outfile.write("PLEASE DOUBLE-CHECK IF THE TEXT MATCHES YOUR EXPECTATIONS")
    outfile.write("TO BE USED FOR EASILY COMMUNICATING A SET-UP PRE-PUBLICATION")
    outfile.write(preamble)
    outfile.write("")
    outfile.write(paragraph)
    outfile.write(text)
    outfile.write("")
    outfile.write("")

with open(bibname, "w") as bibfile:
    bibfile.write(bib)
