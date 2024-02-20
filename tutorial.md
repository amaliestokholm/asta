We have a single star, 16 Cyg A, and we want to know its radius and mass.

### A simple set-up
We have a given project directory. 

I like association-driven project names, so I would personally call it `projectname = 'odette'` after the lead role in Swan Lake - inspired by 16 Cyg A's celestial position in the northern constellation Cygnus (the Swan).

Here I'll describe how I would normally set-up.

I make a folder with my projectname and enter my new project folder
```
PROJECTNAME="odette" && mkdir $PROJECTNAME && cd $PROJECTNAME
```

##### Observational data
In my project folder, I make a folder that contains my observational data.
```
mkdir data && cd data/
```
I make a small divide of my observational constraints: (1) simple observational constraints that is saved in a table and (2) the individual mode frequencies for my target(s) if that is available and something I like to fit.
##### Simple observational constraints
Inside the `data` folder, I put my table (here I'll call it `input.ascii`) that contains a table with data such as the `starid`, effective temperature `teff`, metallicity `MeH`, and other observational data where I have a central value and an uncertainty (`numax`, `dnu`, `LPhot`, etc.). The file could look like this:
```input.ascii
# starid    RA      DEC     numax   enumax  dnu     ednu    teff    eteff   feh     efeh    logg   elogg
16CygA      295.45  50.53   2188.5  20.0    103.3   0.50    5825    50      0.10    0.026   4.33   0.07
```
###### Extra: Individual mode frequencies
To keep the structure simple if I have multiple targets, I prefer to save file(s) containing the individual mode frequencies in a subfolder of `data` that I here will call `freqs`
```
cd data/ && mkdir freqs
```

Suppose you have a text file `16cyga_freqs.ascii` with your individual mode frequencies ala
```16cyga_freqs.ascii
#order, degree, frequency, error (symmetric)
11    1     1334.29   1.01
12    0     1390.81   0.86
12    1     1437.39   0.45
12    2     1487.83    0.71
13    0     1495.05   0.235
13    1     1542.06   0.14
13    2     1590.37   0.39
14    0     1598.69   0.07
14    1     1645.14     0.11
14    2     1693.94   0.19
15    0     1700.95   0.10
15    1     1747.20   0.08
15    2     1795.84   0.13
16    0     1802.35   0.08
```
For the time being, BASTA requires the input of individual mode frequencies to be in an xml file. This is easily generated like so: Activate your `BASTA` environment, navigate to `{$PROJECTNAME}/data/freqs` and run 

```
source ~/BASTA/venv/bin/activate && cd freqs && python -c "import os; from basta import fileio; starid = '16CygA'; fileio.freqs_ascii_to_xml('.', starid, freqsfile=f'{starid}_freqs.ascii')"
```

You do not need to have the radial orders n (`order`) in your file. If your file instead looks like this:
```
# degree, frequency, error
1     1334.29   1.01
0     1390.81   0.86
1     1437.39   0.45
2     1487.83   0.71
0     1495.05   0.24
1     1542.06   0.14
2     1590.37   0.39
0     1598.69   0.07
1     1645.14   0.11
2     1693.94   0.19
0     1700.95   0.10
1     1747.20   0.08
2     1795.84   0.13
0     1802.35   0.08
```
You can still run the same command as before and it will approximate the radial orders for you.



#### Setting up multiple BASTA runs
In the main project folder, I now create a file called `prep_xml.py` in which I set up the grid.
```python
"""
Build cases for BASTA runs using traillib in asta.
"""
import os
import argparse
from asta import traillib

parser = argparse.ArgumentParser()
parser.add_argument("--no-plots", action="store_true")
parser.add_argument("--optional-outputs", action="store_true")
parser.add_argument("action", choices=["xml", "combine", "plot"])


def main() -> None:
    traillib.ignore_astropy_warnings()
    args = parser.parse_args()

    projectname = "odette"
    here = os.path.abspath('.')
    datfile = os.path.join(here, 'data/input.ascii')
    assert os.path.exists(datfile)
    gridlocs = [
            os.path.expanduser('~/BASTA/grids/Garstec_16CygA.hdf5'),
    ]
    assert os.path.isfile(datfile)

    priors = {
        "IMF": "salpeter1955",
        "Teff": {"abstol": "300"},
    }

    # key: value -> colname in .dat file: colname in BASTA.
    inputparams = {
        "starid": "starid",
        "RA": "RA",
        "DEC": "DEC",
        "numax": "numax",
        "numax_err": "enumax",
        "dnu": "dnu",
        "dnu_err": "ednu",
        "Teff": "teff",
        "Teff_err": "eteff",
        "FeH": "feh",
        "FeH_err": "efeh",
        "logg": "logg",
        "logg_err": "elogg",
    }

    outparams = (
        "massfin",
        "radPhot",
        "rho",
        "logg",
        "age",
        "LPhot",
        "Teff",
        "FeH",
        "MeH",
        "dnuSer",
        "numax",
   )

    plotparams = [
        "age",
        "radPhot",
        "massini",
        "Teff",
        "FeH",
        "MeH",
        "LPhot",
        "dnuSer",
        "numax",
        "logg",
    ]

    trail = traillib.init_trail(
        optionaloutputs=args.optional_outputs,
        no_plots=args.no_plots,
        plotparams=plotparams,
        outparams=outparams,
        priors=priors,
        bayweights=True,
        nparts=1,
    )
    
   trail.inputtable = trail.find_inputtable(
        datfile,
        inputparams=inputparams,
        delimiter=",",
        missingval=-9999,
        dustframe="icrs",
    )

    trail.grid = trail.find_grid(
        *gridlocs,
        solarmodel=True,
        solarnumax=3090,
        solardnu=135.1,
        odea=(0.0, 0.0, 0.0, 0.0),
    )

    trail.add_case(
        "01",
        plotname='Teff, FeH',
        fitparams=("Teff", "FeH"),
    )
    trail.add_case(
        "02",
        plotname='Teff, FeH, logg',
        fitparams=("Teff", "FeH", "logg",),
    )
    trail.add_case(
        "03",
        plotname='Teff, FeH, l=0 freqs (BG14)',
        fitparams=("Teff", "FeH", "freqs"),
        freqparams = {
            "freqpath": os.path.join(here, 'data/freqs'),
            "fcor": "BG14",
            "correlations": False,
            "dnufrac": 0.15,
            "onlyradial": True,
            }
    )
    trail.add_case(
        "04",
        plotname='Teff, FeH, freqs (BG14)',
        fitparams=("Teff", "FeH", "freqs"),
        freqparams = {
            "freqpath": os.path.join(here, 'data/freqs'),
            "fcor": "BG14",
            "correlations": False,
            "dnufrac": 0.15,
            "onlyradial": False,
            }
    )

    if args.action == "xml":
        trail.generate_xmls()
    elif args.action == "combine":
       trail.merge_results_together()
    elif args.action == "plot":
        trail.make_comparisonsplot()


if __name__ == "__main__":
    main()
```
This script does all the `BASTA` set-upping for me. Let's go through it from the top.
I'll ignore the lines you also need to ignore:


```python
projectname = "odette"
```
This sets your projectname so everything is named as you like it to be named.

```python3
datfile = os.path.join(here, 'data/input.ascii')
```
This should point to the table file you made earlier in [##### Simple observational constraints].

```
gridlocs = [
    '~/BASTA/grids/Garstec_16CygA.hdf5',
]
```
This should point to the location of the grid you are fitting to. Here we are fitting to the small tutorial grid that ships with `BASTA`.

```python3
inputparams = {
	"starid": "starid",
	"Teff": "Teff",
	"Teff_err": "Teff_err",
	"MeH": "MeH",
	"MeH_err": "MeH_err",
	"LPhot": "LPhot",
	"LPhot_err": "LPhot_err",
	"dnu": "dnu",
	"dnu_err": "dnu_err",
	"numax": "numax",
	"numax_err": "numax_err",
}
```
This mapping shows the mapping between column names in `input.ascii` (keys) and the internal `BASTA` names (values).

```
trail = traillib.init_trail(
        optionaloutputs=args.optional_outputs,
        no_plots=args.no_plots,
        plotparams=plotparams,
        outparams=outparams,
        priors=priors,
        bayweights=True,
        nparts=1,
    )
```
This initialises the run. One important parameter is the `nparts` argument. As we in this case only have 1 star, it does not matter much -- but in cases with many targets, you can split up the total sample into `nparts` smaller parts and make use of the parallelisation to speed up the computation. 

Every `trail.add_case` case is a separate `BASTA` fitting case. Here I am doing 4 fitting cases. In each one I specify:
- A fitting number, e.g. `01` or `03`.
- `plotname`: A single string that in human-readable text explains what the case is.
- `fitparams`: a tuple of all the things I want to fit.
- If `freqs` is among the `fitparams`: `freqparams` is defined. This contains additional fitting settings such as:
	- `freqpath`: the path to the frequency file you defined in [Extra: Individual mode frequencies]
	- `fcor`: which type of surface effect correction I want to use. See the `BASTA` documentation for the possible values.
	- `onlyradial` is a flag that if True makes sure you only fit the l=0 modes. It can be useful to see how the solution might change with the inclusion of the non-radial modes and it is a useful check to see that you can find similar solutions with and without the non-radial modes.

#### How 
Okay, now it is time to run BASTA.
Make sure you have sourced the virtual environment and then run the `prep_xml` script with the argument `xml`.
```
source ~/BASTA/venv/bin/activate && python3 prep_xml xml
```
This should give you an output of the script writing 4 different xml files containing 1 star.

Now let's give them to `BASTA`:
```
cd xmlinputs && BASTAmultirun current --seed 11
```
So we are moving to the newly-created `xmlinputs` folder. In here you can find the xmlfiles in a folder that contains todays date. Furthermore, there is a symlink called `current` pointing to the most recent run. We then run `BASTAmultirun` on the symlink and for good measure we also set the random seed used in `BASTA`, so it is easy to recreate the exact output.

Now wait for this to finish. When it is done, you can go back to your primary folder and have a look at the results in the `results` folder. In here you can have a look at the plots for the single star in each case as well as read the log file.

In order to make it a bit easier for post-processing you can also do the two following steps.
```
source ~/BASTA/venv/bin/activate && python3 prep_xml combine && python3 prep_xml plot 
```
You know have two additional folders: `mainresults`. 

`mainresults` contains subfolders corresponding to each fitting case, i.e. `case01`, `case02` etc. Inside it contains an ascii and a fits file with a summary of the results of that case. This is especially useful if you have many targets and make use of the `nparts` argument as described above.

`plots` now contains an `overview` figure that visually compares the results for each run. This does not tell the full story -- have a look at log files and an inspection of e.g. echelle diagrams  and the individual corner plots for that -- but it makes it easier to spot if one case or a series of cases that might make use of e.g. the same observational constraint are off or different compared to the rest.

Great -- congratulation! You just successfully set up `BASTA` and fitted your target! You can now add more fitting cases to your `prep_xml` script, make new plots using the summarising results files in `mainresults` or whatever you fancy!