"""mplrc -- Set some global matplotlib parameters."""

import matplotlib as mpl


def load_rcparams(style='running') -> None:

    mpl.rcParams['figure.constrained_layout.use'] = True
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.linewidth'] = 0.8
    mpl.rcParams['grid.alpha'] = 0.3

    # Running figures
    if style == 'running':
        pass

    # Report figures
    if style == 'report':
        report_tw = 5.90666  # \printinunitsof{in}\prntlen{\textwidth}
        mpl.rcParams['figure.figsize'] = (report_tw, report_tw/1.618033989)

        mpl.rcParams['figure.dpi'] = 150
        mpl.rcParams['pgf.texsystem'] = "pdflatex"

        mpl.rcParams['mathtext.fontset'] = 'stix'
        mpl.rcParams['font.family'] = 'serif'
        mpl.rcParams['font.serif'] = ['STIX Two Text'] + mpl.rcParams['font.serif']
        mpl.rcParams['font.size'] = 11

        # Those sizes are relative to font.size
        mpl.rcParams['axes.titlesize'] = "medium"
        mpl.rcParams['axes.labelsize'] = "small"
        mpl.rcParams['xtick.labelsize'] = "x-small"
        mpl.rcParams['ytick.labelsize'] = "x-small"

    # Presentation figures
    if style == 'slide':
        mpl.rcParams['mathtext.fontset'] = 'stixsans'
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = ['Noto Sans'] + mpl.rcParams['font.sans-serif']
        mpl.rcParams['font.size'] = 15

        # Those sizes are relative to font.size
        mpl.rcParams['axes.titlesize'] = "medium"
        mpl.rcParams['axes.labelsize'] = "medium"
        mpl.rcParams['xtick.labelsize'] = "small"
        mpl.rcParams['ytick.labelsize'] = "small"
