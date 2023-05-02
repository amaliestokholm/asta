import os
import numpy as np
from cycler import cycler
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from astropy.table import Table, unique

# plt.style.use('/home/amalie/cassiopy/basta_standard.mplstyle')


def get_colors():
    colors = {"nuindiD": "#D55E00", "nuindi": "#0072B2", "hr7322": "#0072B2"}
    tableau10 = [
        "#006BA4",
        "#FF800E",
        "#ABABAB",
        "#595959",
        "#5F9ED1",
        "#C85200",
        "#898989",
        "#A2C8EC",
        "#FFBC79",
        "#CFCFCF",
    ]
    return colors, tableau10


def get_dicts():
    paramsdict = {
        "massfin": "mass",
        "Teff": "teff",
        "dnuSer": "dnu",
        "dnufit": "dnu",
        "numax": "numax",
        "age": "age",
        "FeH": "feh",
        "LPhot": "LPhot",
        "radPhot": "radphot",
        "d02fit": "d02fit",
        "d02mean": "d02mean",
        "R01": "R01",
        "E01": "E01",
    }
    prettypar = {
        "massfin": "$M$\n(M$_{\\odot}$)",
        "radPhot": "$R$\n(R$_{\\odot}$)",
        "FeH": "[Fe/H]\n(dex)",
        "Teff": "$T_{\\mathrm{eff}}$\n(K)",
        "age": "Age\n(Gyr)",
        "Rbcz": "R$_\\mathrm{bcz}$\n(r/R$_\\mathrm{phot}$)",
        "LPhot": "$L$\n(L$_\\odot$)",
        "rho": "$\\rho$\n(g/cm$^3$)",
        "FeHini": "[Fe/H]$_{\\mathrm{ini}}$\n(dex)",
        "yini": "Y$_{\\mathrm{ini}}$",
        "zini": "Z$_{\\mathrm{ini}}$",
        "alphaMLT": "$\\alpha_{\\mathrm{MLT}}$",
        "logg": "$\\log \\, g$\n(dex)",
        "dnuSer": "$\\Delta \\nu_{\\mathrm{Ser}}$\n($\\mu$Hz)",
        "dnufit": "$\\Delta \\nu_{\\mathrm{fit}}$\n($\\mu$Hz)",
        "d02fit": "$\\delta \\nu_{02,\\mathrm{fit}}$\n($\\mu$Hz)",
        "d02mean": "$\\delta \\nu_{02,\\mathrm{mean}}$\n($\\mu$Hz)",
        "R01": "r(01)",
        "E01": "$\delta \epsilon_{01}$",
        "numax": "$\\nu_{\\mathrm{max}}$\n($\\mu$Hz)",
        "parallax": "$\\varpi$",
    }

    return paramsdict, prettypar


def get_name(bastaoutdirs):
    # Make yaxis
    labels = []
    for i in bastaoutdirs:
        name = r""
        if "Teff" in i or "TEFF" in i:
            name += "$T_{\mathrm{eff}}$, "
        if "FeH" in i or "FEH" in i:
            name += "[Fe/H], "
        if "MeH" in i:
            name += "[M/H], "
        if "dnufit" in i or 'DNU' in i:
            name += r"$\Delta \nu_{\mathrm{fit}}$, "
        if "dnuSer" in i or "SER" in i:
            name += r"$\Delta \nu_{\mathrm{Ser}}$, "
        if "numax" in i or "NUM" in i:
            name += r"$\nu_{\mathrm{max}}$, "
        if "LPhot" in i or 'LPHOT' in i:
            name += r"$L$, "
        if "PHASE" in i:
            name += "phase, "
        if "E01" in i or 'EPSDIFF' in i:
            name += r"$\delta \epsilon_{01}$., "
        if "R01" in i or 'RATIO' in i:
            name += r"ratios (01)., "
        if "freqs" in i or 'FREQ' in i:
            name += r"freqs., "
            if 'CUBIC' in i:
                name+='BG14cubic'
            if "minus" in i:
                name = name[:-2]
                name += " - $(l, n) =$"
            if "l0n11" in i:
                name += "($0,11$), "
            if "l1n7n8" in i:
                name += "($1,7$), ($1,8$), "
            if "l1n16" in i:
                name += "($1,16$), "
            if "l2n7" in i:
                name += "($2,7$), "
        if "parallax" in i or "MAG" in i or "TYCHO" in i or "GAIA" in i or "BPRP" in i:
            name += r"$\varpi$+("
            if "Mh_2MASS" in i or "HMAG" in i:
                name += r"$H$, "
            if "Mj_2MASS" in i or "JMAG" in i:
                name += r"$J$, "
            if "Mk_2MASS" in i or "KMAG" in i:
                name += r"$K$, "
            if "Mv_JC" in i:
                name += r"$V$, "
            if "Mv_TYCHO" in i or "VTMAG" in i or "TYCHO" in i:
                name += r"$V_T$, "
            if "Mb_TYCHO" in i or "BTMAG" in i or "TYCHO" in i:
                name += r"$B_T$, "
            if "G_GAIA" in i or "GMAG" in i or "GAIA" in i:
                name += r"$G$, "
            if "RP_GAIA" in i or "RPMAG" in i or "GAIA" in i or "BPRP" in i:
                name += r"$RP$, "
            if "BP_GAIA" in i or "BPMAG" in i or "GAIA" in i or "BPRP" in i:
                name += r"$BP$, "
            name = name[:-2] + ")"
        else:
            name = name[:-2]
        labels.append(name)
    return labels


def compare_different_runs(
    inptable,
    bastaoutdirs,
    idstr,
    parameter_for_comparison,
    mode="",
    sourceid="source_id",
    style=None,
    prior=None,
    includeprior=True,
    plotdir='.',
):
    if style is None:
        prettifys = [True, False]
    elif style == "pretty":
        prettifys = [True]
    else:
        prettifys = [False]

    paramsdict, prettypar = get_dicts()
    colors, tableau10 = get_colors()
    inp = Table.read(inptable, format="ascii.commented_header", delimiter=',')

    labels = get_name(bastaoutdirs)

    stars = inp[sourceid]
    for i, starid in enumerate(stars):
        fname = f"{plotdir}/{str(idstr)}{str(mode)}_{str(starid)}.pdf"
        if str(idstr) in colors.keys():
            c = [colors[str(idstr)]] * len(bastaoutdirs)
            obsc = tableau10[6]
        else:
            c = [tableau10[i % len(tableau10)]] * len(bastaoutdirs)
            obsc = tableau10[i % len(tableau10)]

        with PdfPages(fname) as pdf:
            for prettify in prettifys:
                fig, axs = plt.subplots(
                    1, len(parameter_for_comparison), sharey=True
                )  # , figsize=(4, 3))
                fig.subplots_adjust(wspace=0)
                axs[0].set_yticks(np.arange(len(bastaoutdirs)))
                axs[0].invert_yaxis()
                if prettify:
                    axs[0].set_yticklabels(labels, fontsize="small")
                else:
                    axs[0].set_yticklabels(bastaoutdirs, fontsize="small")

                for p, par in enumerate(parameter_for_comparison):
                    inppar = paramsdict[par]
                    if inppar in inp.columns:
                        obs = inp[inppar][i]
                        if inppar + "_e" in inp.columns:
                            obs_err = inp[inppar + "_e"][i]
                        else:
                            obs_err = inp[inppar + "_err"][i]

                    else:
                        obs = None

                    if not obs is None:
                        axs[p].axvline(x=obs, c=obsc, alpha=0.7, zorder=0)
                        axs[p].axvspan(
                            obs - obs_err,
                            obs + obs_err,
                            color=obsc,
                            alpha=0.3,
                            zorder=0,
                        )
                    xs = []
                    errs = []
                    for j, bastaoutdir in enumerate(bastaoutdirs):
                        if os.path.isdir(bastaoutdir):
                            for f in os.listdir(bastaoutdir):
                                if f.endswith(".ascii") and not f.endswith(
                                    "distance.ascii"
                                ):
                                    asciifile = os.path.join(bastaoutdir, f)

                            try:
                                bastatable = Table.read(
                                    asciifile, format="ascii.commented_header",
                                )
                            except Exception:
                                print("Could not find ascii table in")
                                print(f)
                                raise SystemExit
                            if starid not in bastatable["starid"]:
                                print(f"{starid} not in bastatable")
                                continue
                            else:
                                star = bastatable["starid"] == starid
                            try:
                                x = bastatable[par][star].value[0]
                            except Exception:
                                print("exception")
                                if par == "alphaMLT":
                                    x = bastatable['alpha'][star].value[0]
                                else:
                                    raise
                            if par == "age":
                                x /= 1e3
                            xs.append(x)

                            try:
                                merr = bastatable[par + "_errm"][star].value[0]
                                perr = bastatable[par + "_errp"][star].value[0]
                            except Exception:
                                if par == "alphaMLT":
                                    par = "alpha"
                                    merr = bastatable[par + "_errm"][star].value[0]
                                    perr = bastatable[par + "_errp"][star].value[0]
                                    par = "alphaMLT"
                                else:
                                    raise
                            xerr = np.asarray([merr, perr]).reshape((2,1))
                            if par == "age":
                                xerr[0] /= 1e3
                                xerr[1] /= 1e3
                            errs.extend(xerr)
                            if np.isnan(x):
                                axs[p].errorbar(x, j, marker="x", c="k")
                            else:
                                axs[p].errorbar(
                                    x,
                                    j,
                                    xerr=xerr,
                                    marker=".",
                                    mfc=c[j],
                                    markeredgecolor="k",
                                    fillstyle="full",
                                    markeredgewidth=0.5,
                                    ecolor="k",
                                    elinewidth=2,
                                )
                                axs[p].tick_params(axis="both", labelsize="small")

                    if obs is not None:
                        xs.append(obs)

                    assert len(xs) > 0, par
                    axs[p].set_xlabel(prettypar[par], fontsize="small")
                    if not np.isnan(xs).all():
                        l = np.nanmin(xs) - 1.2 * np.nanmax(errs)
                        u = np.nanmax(xs) + 1.2 * np.nanmax(errs)
                        if not np.isfinite(l) or not np.isfinite(u):
                            axs[p].set_xlim([-1, 1])
                        elif l == u:
                            axs[p].set_xlim([np.nanmin(xs) * 0.9, np.nanmax(xs) * 1.1])
                        else:
                            axs[p].set_xlim([l, u])
                    plt.setp(
                        axs[p].get_xticklabels(),
                        rotation=30,
                        horizontalalignment="right",
                    )
                    fig.align_xlabels(axs[:])

                plt.title(str(starid), fontsize="small")
                pdf.savefig(bbox_inches="tight")
                plt.close()


def plot_comparison_of_different_stars(asciifiles, idstrs, fname):
    tables = [Table.read(asciifile, format="ascii") for asciifile in asciifiles]
    with PdfPages(fname) as pdf:
        f, axs = plt.subplots(1, 2, sharey=True, sharex=True)
        f.subplots_adjust(wspace=0.1)
        tables[0] = unique(tables[0], keys="starid")
        tables[1] = unique(tables[1], keys="starid")
        starids = tables[0]["starid"].data.astype(str)
        axs[0].set_yticks(np.arange(len(starids)))
        axs[0].invert_yaxis()
        starids = [starid.replace("b'", "").replace("'", "") for starid in starids]
        axs[0].set_yticklabels(starids, fontsize="small")
        overlap = np.in1d(tables[0]["starid"].data, tables[1]["starid"].data)
        nooverlap = np.where(~overlap)[0]
        print("Overlap", nooverlap)
        for (table, ax, idstr) in zip(tables, [axs[0], axs[1]], idstrs):
            # Replace nan with something else
            mask = np.isnan(table["massfin"])
            table["massfin"][mask] = 0
            table["massfin_errm"][mask] = 0
            table["massfin_errp"][mask] = 0
            merr = table["massfin_errm"].data
            perr = table["massfin_errp"].data
            massfin = table["massfin"].data
            # Fill in with x's where stars in table 0 are not in table 1
            if ax is axs[1]:
                for extrai in nooverlap:
                    massfin = np.insert(massfin, extrai, 0)
                    merr = np.insert(merr, extrai, 0)
                    perr = np.insert(perr, extrai, 0)
            for i, M in enumerate(massfin):
                xerr = [[merr[i]], [perr[i]]]
                if M == 0:
                    ax.errorbar(M, i, marker="x", c="k")
                else:
                    ax.errorbar(
                        M,
                        i,
                        xerr=xerr,
                        marker="o",
                        markeredgecolor="k",
                        fillstyle="full",
                        markeredgewidth=0.5,
                        ecolor="k",
                        elinewidth=2,
                    )
            ax.text(max(massfin), 0.1, idstr, ha="left", size="small")
            # ax.set_xlim([-0.1, max(massfin) + 0.5])
            ax.set_xlim([0.45, 2.05])
            ax.tick_params(axis="both", labelsize="small")
        f.text(0.5, 0.04, r"Stellar mass (M$_{\odot}$)", ha="center")
        pdf.savefig(bbox_inches="tight")


def plot_comparison_of_different_bastaoutput(
    starids,
    mainfolder,
    dirs,
    idstr,
    parameter_for_comparison,
    genname,
    filename,
    observed,
    observed_err,
    style=None,
    prior=None,
    includeprior=True,
    horizontal=True,
):
    if style is None:
        prettifys = [True, False]
    elif style == "pretty":
        prettifys = [True]
    else:
        prettifys = [False]
    colors, tableau10 = get_colors()
    paramsdict, prettypar = get_dicts()
    bastaoutputs = []
    if prior is not None:
        for bastaoutput in dirs:
            if includeprior:
                if prior in bastaoutput:
                    bastaoutputs.append(bastaoutput)
            else:
                if prior not in bastaoutput:
                    bastaoutputs.append(bastaoutput)

    else:
        bastaoutputs = dirs
    assert len(bastaoutputs) > 0

    labels = get_name(bastaoutdirs)

    for i, starid in enumerate(starids):
        if prior is not None:
            fname = genname + idstr + "_" + str(starid) + "_" + str(prior) + ".pdf"
        else:
            fname = genname + idstr + "_" + str(starid) + ".pdf"
        if idstr == "compare":
            c = [colors[str(name).split("_")[0]] for name in bastaoutputs]
            obsc = tableau10[6]
        else:
            if str(idstr) in colors.keys():
                c = [colors[str(idstr)]] * len(bastaoutputs)
                obsc = tableau10[6]
            else:
                c = [tableau10[i % len(tableau10)]] * len(bastaoutputs)
                obsc = tableau10[i % len(tableau10)]
        if horizontal:
            with PdfPages(fname) as pdf:
                for prettify in prettifys:
                    fig, axs = plt.subplots(
                        1, len(parameter_for_comparison), sharey=True
                    )  # , figsize=(4, 3))
                    fig.subplots_adjust(wspace=0)
                    axs[0].set_yticks(np.arange(len(dirs)))
                    axs[0].invert_yaxis()
                    if prettify:
                        axs[0].set_yticklabels(labels, fontsize="small")
                    else:
                        axs[0].set_yticklabels(bastaoutputs, fontsize="small")
                    for i, par in enumerate(parameter_for_comparison):
                        if not observed[i] is None:
                            axs[i].axvline(x=observed[i], c=obsc, alpha=0.7, zorder=0)
                            axs[i].axvspan(
                                observed[i] - observed_err[i],
                                observed[i] + observed_err[i],
                                color=obsc,
                                alpha=0.3,
                                zorder=0,
                            )
                        xs = []
                        errs = []
                        for j, bastaoutput in enumerate(bastaoutputs):
                            f = os.path.join(mainfolder, bastaoutput)
                            if filename == "all":
                                for context in os.listdir(f):
                                    if content.endswith(".ascii"):
                                        print("uses the file ending in ascii")
                                        f = os.path.join(f, content)
                            else:
                                f = os.path.join(f, filename)

                            try:
                                bastatable = Table.read(f, format="ascii")
                            except Exception:
                                print("Could not find ascii table in")
                                print(f)
                                raise SystemExit
                            star = bastatable["starid"] == starid
                            try:
                                x = bastatable[par][star]
                            except Exception:
                                print(bastaoutput)
                                if par == "alphaMLT":
                                    par = "alpha"
                                    x = bastatable[par][star]
                                    par = "alphaMLT"
                                else:
                                    raise
                            if par == "age":
                                x /= 1e3
                            xs.append(x)
                            try:
                                merr = bastatable[par + "_errm"][star]
                                perr = bastatable[par + "_errp"][star]
                            except Exception:
                                if par == "alphaMLT":
                                    par = "alpha"
                                    merr = bastatable[par + "_errm"][star]
                                    perr = bastatable[par + "_errp"][star]
                                    par = "alphaMLT"
                                else:
                                    raise
                            xerr = [merr, perr]
                            if par == "age":
                                xerr[0] /= 1e3
                                xerr[1] /= 1e3
                            errs.append(np.nanmin(xerr))
                            if np.isnan(x[0]):
                                axs[i].errorbar(x, j, marker="x", c="k")
                            else:
                                axs[i].errorbar(
                                    x,
                                    j,
                                    xerr=xerr,
                                    marker=".",
                                    mfc=c[j],
                                    markeredgecolor="k",
                                    fillstyle="full",
                                    markeredgewidth=0.5,
                                    ecolor="k",
                                    elinewidth=2,
                                )
                                axs[i].tick_params(axis="both", labelsize="small")
                        axs[i].set_xlabel(prettypar[par], fontsize="small")
                        axs[i].set_xlim(
                            [np.nanmin(xs) - 1.2 * np.nanmax(errs), np.nanmax(xs) + 1.2 * np.nanmax(errs)]
                        )
                        plt.setp(
                            axs[i].get_xticklabels(),
                            rotation=30,
                            horizontalalignment="right",
                        )
                    fig.align_xlabels(axs[:])
                    pdf.savefig(bbox_inches="tight")

        else:
            with PdfPages(fname) as pdf:
                for prettify in prettifys:
                    fig, axs = plt.subplots(
                        len(parameter_for_comparison), 1, sharex=True
                    )  # , figsize=(2, 4))
                    fig.subplots_adjust(hspace=0)
                    axs[0].set_xticks(np.arange(len(dirs)))
                    axs[0].invert_xaxis()
                    if prettify:
                        axs[0].set_xticklabels(labels, fontsize="small")
                    else:
                        axs[0].set_xticklabels(bastaoutputs, fontsize="small")
                    for i, par in enumerate(parameter_for_comparison):
                        if not observed[i] is None:
                            axs[i].axhline(y=observed[i], c=obsc, alpha=0.7, zorder=0)
                            axs[i].axhspan(
                                observed[i] - observed_err[i],
                                observed[i] + observed_err[i],
                                color=obsc,
                                alpha=0.3,
                                zorder=0,
                            )
                        xs = []
                        errs = []
                        for j, bastaoutput in enumerate(bastaoutputs):
                            f = os.path.join(mainfolder, bastaoutput)
                            if filename == "all":
                                for context in os.listdir(f):
                                    if content.endswith(".ascii"):
                                        print("Uses the file ending in ascii")
                                        f = os.path.join(f, content)
                            else:
                                f = os.path.join(f, filename)
                            try:
                                bastatable = Table.read(f, format="ascii")
                            except Exception:
                                print("Could not find ascii table in")
                                print(f)
                                raise SystemExit
                            star = bastatable["starid"] == starid
                            try:
                                x = bastatable[par][star]
                            except Exception:
                                print(bastaoutput)
                                if par == "alphaMLT":
                                    par = "alpha"
                                    x = bastatable[par][star]
                                    par = "alphaMLT"
                                else:
                                    raise
                            if par == "age":
                                x /= 1e3
                            xs.append(x)
                            try:
                                merr = bastatable[par + "_errm"][star]
                                perr = bastatable[par + "_errp"][star]
                            except Exception:
                                if par == "alphaMLT":
                                    par = "alpha"
                                    merr = bastatable[par + "_errm"][star]
                                    perr = bastatable[par + "_errp"][star]
                                    par = "alphaMLT"
                                else:
                                    raise
                            xerr = [merr, perr]
                            if par == "age":
                                xerr[0] /= 1e3
                                xerr[1] /= 1e3
                            errs.append(np.nanmin(xerr))
                            if np.isnan(x[0]):
                                axs[i].errorbar(j, x, marker="x", c="k")
                            else:
                                axs[i].errorbar(
                                    j,
                                    x,
                                    yerr=xerr,
                                    marker=".",
                                    mfc=c[j],
                                    markeredgecolor="k",
                                    fillstyle="full",
                                    markeredgewidth=0.5,
                                    ecolor="k",
                                    elinewidth=2,
                                )
                        # axs[i].tick_params(axis='both', labelsize='small')
                        axs[i].set_ylabel(prettypar[par], fontsize="small")
                        axs[i].set_ylim(
                            [np.nanmin(xs) - 1.2 * np.nanmax(errs), np.nanmax(xs) + 1.2 * np.nanmax(errs)]
                        )
                        plt.setp(
                            axs[i].get_xticklabels(),
                            rotation=30,
                            horizontalalignment="right",
                            fontsize="small",
                        )
                    fig.align_xlabels(axs[:])
                    pdf.savefig(bbox_inches="tight")
