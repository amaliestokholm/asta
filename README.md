# asta

This repository contains helper functions for running [BASTA](https://github.com/BASTAcode/BASTA). 

It can be useful (i) if you run BASTA for many stars at a time using BASTAmultirun or (ii) if you run many different fitting cases or other set-ups for a sample of stars.

You can get this locally using pip
```
pip install https://github.com/amaliestokholm/asta/archive/main.zip
```
You would also need numpy, astropy, and matplotlib.

You would then have access to the newest version of `traillib.py`.
In this script, you can find the commands for creating xml files easily using `init_trail`, `add_case` and `generate_xmls()`.


## Quick guide: How to use it
For a full guide, see [tutorial.md](tutorial.md)

You set-up the `prep_xml.py` as you prefer (see tutorial or https://cloud.phys.au.dk/nextcloud/index.php/s/gfq2RsbDnDy8a2E or https://cloud.phys.au.dk/nextcloud/index.php/s/Z7tTQFDE3m7e3Mc for inspiration).

Then run
```
python3 prep_xml.py xml
```
to generate the xmls.

When you are ready to run your fits, you do the following from your project directory.
```
source ../../../BASTA/venv/bin/activate 
cd ./xmlinput
BASTAmultirun current
```
Notice how simple the last line is. Current is a symlink set-up in `traillib` so you do not have to worry about the directory names.

After the run is complete, you can run
```
python3 prep_xml.py combine
```
To combine your resultsfiles and get an overview of how many stars resulted in `nan` outputs.
The stars with `nan` output will then be written to a new xml file, so if they were caused by e.g. server downtime, they can easily be fitted anew.

For visual inspection, make plots using
```
python3 prep_xml.py plot
```

### Many stars? Divide your sample into smaller chunks!
Depending on your value of `nparts` in `trail = traillib.init_trail`, it divides your sample into smaller chunks for easier parallelization. Don't worry - the script can merge your resultsfiles afterwards so you only have 1 file for each fitting case in the end.

### Run a test with a few stars!
If you want to just do a test run with a few stars, I suggest you add the line
```python3
trail.quicktest_onlystars = 3
```
to your `prep_xml.py`.

### Use different reddening values
Want to test the effect of including the E(B-V) from the dustmap?
You can add your own column of E(B-V) values to your input `.ascii` file. This could be named 'EBV'.

In your `prep_xml`.py, you need to not include `"EBV": "EBV",` in your `inputparams`, but instead when you set-up your cases, you can do as follows:
```python3
    trail.add_case(
        "01",
        fitparams=("Teff", "MeH", "parallax"),
        extraplotparams=["distance"],
        filters=["RP_GAIA"],
        ignoreebv=True,
        )
    trail.add_case(
        "02",
        fitparams=("Teff", "MeH", "parallax",),
        extraplotparams=["distance"],
        filters=["RP_GAIA"],
        ignoreebv=False,
        )
```
