import datetime
import glob
import os
import re
import string
from typing import Any, Sequence

import numpy as np
from astropy.table import Table, join, vstack


def ignore_astropy_warnings() -> None:
    import warnings

    from astropy.utils.exceptions import AstropyWarning

    warnings.filterwarnings("ignore", category=AstropyWarning)


Odea = tuple[float, float, float, float]


class Grid:
    filename: str

    solarmodel: bool
    solarnumax: float
    solardnu: float
    odea: Odea


class Inputtable:
    filename: str

    inputparams: dict[str, str]
    delimiter: str
    missingval: float
    dustframe: str


class Case:
    casename: str
    inputtable: Inputtable
    inputparams: dict[str, str]
    plotparams: list[str]
    fitparams: list[str]
    filters: list[str]


class Trail:
    optionaloutputs: bool
    no_plots: bool
    plotparams: list[str]
    outparams: list[str]
    priors: dict[str, Any]
    bayweights: bool

    inputtable: Inputtable | None = None
    grid: Grid | None = None
    cases: list[Case]
    nparts: int
    runid: str
    xml_file_count: int

    quicktest_onlystars: int | None = None

    def find_grid(
        self,
        *filenames: str,
        solarmodel: bool,
        solarnumax: float,
        solardnu: float,
        odea: Odea,
    ) -> Grid:
        assert filenames
        filename = find_file(filenames)
        if filename is None:
            raise SystemExit("Couldn't find grid at %r" % (filenames[0],))
        grid = Grid()
        grid.filename = filename
        grid.solarmodel = solarmodel
        grid.solarnumax = solarnumax
        grid.solardnu = solardnu
        grid.odea = odea
        return grid

    def find_inputtable(
        self,
        *filenames: str,
        inputparams: dict[str, str],
        delimiter: str,
        missingval: float,
        dustframe: str,
    ) -> Inputtable:
        assert filenames
        filename = find_file(filenames)
        if filename is None:
            raise SystemExit("Couldn't find input table at %r" % (filenames[0],))
        inputtable = Inputtable()
        inputtable.filename = filename
        inputtable.inputparams = {**inputparams}
        inputtable.delimiter = delimiter
        inputtable.missingval = missingval
        inputtable.dustframe = dustframe
        return inputtable

    def add_case(
        self,
        casename: str,
        *,
        inputtable: Inputtable | None = None,
        inputparams: dict[str, str] | None = None,
        extraplotparams: Sequence[str] | None = None,
        fitparams: Sequence[str],
        filters: Sequence[str] | None = None,
    ) -> None:
        if not isinstance(casename, str) or not re.match(r"^[0-9][0-9]$", casename):
            raise Exception("Case name must be a string of two digits")

        if inputtable is None:
            assert self.inputtable is not None
            inputtable = self.inputtable

        if inputparams is None:
            inputparams = {}
        inputparams = {**inputtable.inputparams, **inputparams}

        if extraplotparams is None:
            extraplotparams = ()

        if filters is None:
            filters = ()

        thecase = Case()
        thecase.casename = casename
        thecase.inputtable = inputtable
        thecase.inputparams = inputparams
        thecase.plotparams = [*self.plotparams, *extraplotparams]
        thecase.fitparams = list(fitparams)
        thecase.filters = list(filters)
        self.cases.append(thecase)

    def generate_xmls(self) -> None:
        mkdir_if_not_exists("xmlinputs")
        mkdir_if_not_exists("results")
        self._choose_available_runid()
        remove_if_exists("xmlinputs/current")
        os.symlink(self.runid, "xmlinputs/current")
        self.xml_file_count = 0
        for thecase in self.cases:
            self._generate_xml(thecase)
        if not self.xml_file_count:
            raise SystemExit("All stars already done")

    def _choose_available_runid(self) -> None:
        today = datetime.datetime.now().strftime("%Y%m%d")
        for s in string.ascii_lowercase:
            dirname = f"xmlinputs/{today}{s}"
            if os.path.exists(dirname):
                continue
            self.xmlinput_dir = dirname
            self.runid = os.path.basename(dirname)
            break
        else:
            raise Exception("Failed to create xmlinput dir")
        self.results_dir = f"results/{self.runid}"

    def _generate_xml(self, thecase: Case) -> None:
        from basta.xml_create import generate_xml

        inputtable = thecase.inputtable
        grid = self.grid
        assert grid is not None

        oldresults = sorted(glob.glob(f"results/*/*/case{thecase.casename}*.ascii"))
        donestars: list[str] = []
        for oldresultfile in oldresults:
            if os.stat(oldresultfile).st_size == 0:
                continue
            with open(oldresultfile) as fp:
                header = fp.readline().lstrip(" #").split()
                assert header[0] == "starid", header
                assert len(header) > 1
                for line in fp:
                    starid, *values = line.split()
                    assert values
                    if values[0] == "nan":
                        continue
                    # check that it's a float value
                    float(values[0])
                    donestars.append(starid)
        donestars_set = set(donestars)

        inp = Table.read(inputtable.filename, format="ascii.commented_header")
        if self.quicktest_onlystars is not None:
            inp = inp[: self.quicktest_onlystars]
        first_column_name = list(inp.columns)[0]
        first_column = inp.columns[first_column_name]
        assert first_column.name == thecase.inputparams["starid"]
        donemask = np.array([starid in donestars_set for starid in first_column])

        inp = inp[~donemask]
        if len(inp) == 0:
            print(f"Case {thecase.casename}: All stars already done")
            return
        nparts = min(len(inp), self.nparts)
        tables = [
            inp[len(inp) * i // nparts : len(inp) * (i + 1) // nparts]
            for i in range(nparts)
        ]
        assert sum(len(t) for t in tables) == len(inp)

        reverse_inputcols: dict[str, str] = {}
        for bastacol, inputcol in thecase.inputparams.items():
            if inputcol in reverse_inputcols:
                raise Exception(
                    "Input column %r is mapped to both %r and %r - this is not supported"
                    % (inputcol, reverse_inputcols[inputcol], bastacol)
                )
            reverse_inputcols[inputcol] = bastacol
        assert len(reverse_inputcols) == len(thecase.inputparams)

        inputparams: list[str] = []
        for colname in inp.columns:
            try:
                bastacol = reverse_inputcols[colname]
            except KeyError:
                assert colname not in thecase.inputparams
                bastacol = colname
            inputparams.append(bastacol)
        missing = thecase.inputparams.keys() - set(inputparams)
        if missing:
            raise Exception(
                "Missing columns in input: %s"
                % (
                    ", ".join(
                        "%s (%s)" % (thecase.inputparams[k], k)
                        for k in thecase.inputparams
                        if k in missing
                    )
                )
            )
        for i, tab in enumerate(tables):
            basename = f"case{thecase.casename}_part{i+1}"
            asciifile = f"{self.xmlinput_dir}/{basename}.dat"
            mkdir_if_not_exists(self.xmlinput_dir)
            mkdir_if_not_exists(self.results_dir)
            outputpath = os.path.join(self.results_dir, f"case{thecase.casename}/")
            mkdir_if_not_exists(outputpath)
            delimiter = ","
            tab.write(asciifile, format="ascii.commented_header", delimiter=delimiter)
            missingval = inputtable.missingval
            outputfile = f"{basename}.ascii"
            xmlfile = f"{self.xmlinput_dir}/{basename}.xml"

            print(f"Write {xmlfile} containing {len(tab)} stars")
            xml = generate_xml(
                gridfile=os.path.abspath(grid.filename),
                asciifile=os.path.abspath(asciifile),
                outputpath=os.path.abspath(outputpath),
                odea=grid.odea,
                params=inputparams,
                fitparams=thecase.fitparams,
                outparams=self.outparams,
                missingval=missingval,
                priors=self.priors,
                bayweights=self.bayweights,
                optionaloutputs=self.optionaloutputs,
                outputfile=outputfile,
                kielplots=not self.no_plots,
                cornerplots=[] if self.no_plots else self.plotparams,
                dustframe=inputtable.dustframe,
                filters=thecase.filters,
                solarmodel=True,
                sunnumax=grid.solarnumax,
                sundnu=grid.solardnu,
                delimiter=delimiter,
            )

            self.xml_file_count += 1
            with open(xmlfile, "w") as inpfile:
                print(xml, file=inpfile)

    # Merge the results ascii files together
    def merge_results_together(self, verbose=False):
        for thecase in self.cases:
            base = f"mainresults/case{thecase.casename}/case{thecase.casename}"
            mkdir_if_not_exists("mainresults")
            mkdir_if_not_exists(f"mainresults/case{thecase.casename}")
            merge_results_together(
                inputs=sorted(glob.glob(f"results/*/*/case{thecase.casename}*.ascii")),
                outputs=[f"{base}.ascii", f"{base}.fits"],
                verbose=verbose,
            )


def init_trail(
    *,
    optionaloutputs: bool,
    no_plots: bool,
    plotparams: Sequence[str] = (),
    outparams: Sequence[str],
    priors: dict[str, Any],
    bayweights: bool,
    nparts: int | None = None,
) -> Trail:
    if nparts is None:
        nparts = 1
    assert isinstance(nparts, int)
    assert nparts >= 1
    trail = Trail()
    trail.optionaloutputs = optionaloutputs
    trail.no_plots = no_plots
    trail.plotparams = list(plotparams)
    trail.outparams = list(outparams)
    trail.priors = {**priors}
    trail.bayweights = bayweights
    trail.cases = []
    trail.nparts = nparts
    return trail


def find_file(filenames: Sequence[str]) -> str | None:
    for f in filenames:
        if os.path.isfile(f):
            return f

    # Not found - but if one of them is e.g. a directory, raise an Exception
    for f in filenames:
        if os.path.exists(f):
            raise Exception("%r is not a regular file" % (f,))

    # Not found
    return None


def mkdir_if_not_exists(f: str) -> None:
    try:
        os.mkdir(f)
    except FileExistsError:
        pass


def remove_if_exists(f: str) -> None:
    try:
        os.remove(f)
    except FileNotFoundError:
        pass


def merge_results_together(
    *, inputs: list[str], outputs: list[str], verbose: bool = False
) -> None:
    if not inputs:
        raise Exception("At least one input file must be specified")
    if not outputs:
        raise Exception("At least one output file must be specified")
    for o in outputs:
        if not o.endswith((".fits", ".ascii")):
            raise Exception("Output must end in .fits or .ascii: %r" % (o,))

    resultfiles = [f for f in inputs if os.stat(f).st_size > 0]
    if not resultfiles:
        raise Exception("All input files are empty")
    resulttables: list[Table] = []
    for rf in resultfiles:
        try:
            resulttables.append(Table.read(rf, format="ascii"))
        except Exception:
            print(rf)
            raise
    a = resulttables[0]
    verbose and print(
        "First read %s: %s stars, %s with nans"
        % (resultfiles[0], len(a), np.sum(~np.isfinite(a.columns[1])))
    )
    for i, b in enumerate(resulttables):
        if i == 0:
            continue

        b = b[np.argsort(b["starid"])]

        # Find indexes of rows in a with nan values
        nanindex = (~np.isfinite(a.columns[1])).nonzero()[0]
        isin = np.isin(a["starid"][nanindex], b["starid"])
        # verbose and print(a["starid"][nanindex[~isin]])
        nanindex_in_b = nanindex[isin]
        nanids = a["starid"][nanindex_in_b]
        nanid_in_b = np.searchsorted(b["starid"], nanids)
        a[nanindex_in_b] = b[nanid_in_b]

        newstars = b[~np.isin(b["starid"], a["starid"])]
        verbose and print(
            "Now read %s: %s stars, %s with nans, overwrite %s, append %s"
            % (
                resultfiles[i],
                len(b),
                np.sum(~np.isfinite(b.columns[1])),
                isin.sum(),
                len(newstars),
            )
        )
        if len(newstars) > 0:
            a = vstack([a, newstars], join_type="exact", metadata_conflicts="error")

    stillbadstars = ~np.isfinite(a["massfin"])
    stillbadstarids = a["starid"][stillbadstars]
    print(
        "%s: %s with nans out of %s: %s"
        % (
            outputs[0],
            len(stillbadstarids),
            len(a),
            ", ".join(map(str, stillbadstarids.data[:10])),
        )
    )
    verbose and print(stillbadstarids)
    # a = join(inptable, a, keys=["starid"], join_type="left")
    for o in outputs:
        if o.endswith(".fits"):
            a.write(o, format="fits", overwrite=True)
        elif o.endswith(".ascii"):
            a.write(o, format="ascii.commented_header", overwrite=True)
        else:
            raise Exception("unhandled extension: %r" % (o,))
