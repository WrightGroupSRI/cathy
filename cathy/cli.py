"""Command Line Interface for cathy."""
import collections
import datetime
import functools
import logging
import os
import random
import re
from functools import reduce
from operator import getitem
from pathlib import Path
from pprint import pprint

import pydicom as dicom
import click
import click_log
import enlighten
import numpy
import pandas as pd
from matplotlib import animation, pyplot
from datetime import datetime
import os

import catheter_ukf
import catheter_utils.cathcoords
import catheter_utils.geometry
import catheter_utils.localization
import catheter_utils.projections
import catheter_utils.metrics
import dicom_art
from get_gt import __main__ as get_gt
from . import art

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# set up logging

click_log.basic_config()
logger = logging.getLogger(__name__)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# helpers


class ExistingDirectory(click.Path):
    """A click.Path with customized default args."""

    def __init__(self, *args, **kwargs):
        kwargs.update(file_okay=False, dir_okay=True, exists=True)
        super().__init__(*args, **kwargs)


class Ints(click.ParamType):
    """Options strings will be an equals sign followed by one or more
    non-negative integers, separated by commas and dashes. Dashes indicate
    ranges (e.g., "1-3" -> [1, 2, 3]), and commas
    separate values and ranges.
    """

    name = "INTs"

    def convert(self, value, param, ctx):
        if value is None or len(value) == 0:
            return None
        if value[0] == '=':
            value = value[1:-1]

        range_re = re.compile(r"(\d+)-(\d+)")
        values = []
        for part in value.split(","):
            try:
                values.append(int(part))
            except ValueError:
                range_match = range_re.match(part)
                if range_match:
                    values.extend(range(int(range_match.group(1)), 1 + int(range_match.group(2))))
                else:
                    self.fail("unsupported option format: {}".format(value), param, ctx)

        return sorted(values)


class ChoiceOrInts(Ints):
    def __init__(self, choices):
        self.choices = set([choice.lower() for choice in choices])

    @property
    def name(self):
        return "[{}|{}]".format("|".join(self.choices), super().name)

    def convert(self, value, param, ctx):
        lvalue = value.lower()
        if lvalue in self.choices:
            return lvalue
        else:
            return super().convert(value, param, ctx)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cathy root


@click.group()
def cathy():
    """Command line tools for processing catheter data."""
    pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cathy apply_ukf


@cathy.command()
@click.argument("src", type=ExistingDirectory())
@click.argument("dst", type=ExistingDirectory())
@click.option("--distal", "-d", default=4, help="Distal coil index.")
@click.option("--proximal", "-p", default=5, help="Proximal coil index.")
@click_log.simple_verbosity_option()
def apply_ukf(src, dst, distal, proximal):
    """Create UKF smoothed versions of cathcoords files."""

    # - - - - - - - - -
    # 1. Discover data

    # Find cathcoords files in the input and output directories.
    src_info = catheter_utils.cathcoords.discover_files(src)
    dst_info = catheter_utils.cathcoords.discover_files(dst)

    # Is there any data?
    if len(src_info) == 0:
        logger.info("No cathcoords found in input. Done!")
        return

    logger.info("Found {} cathcoords recordings".format(len(src_info)))

    # Are the requested coils found in the available data?
    coils = set.intersection(*[set(r.keys()) for r in src_info.values()])
    requested = set((distal, proximal))
    notfound = requested.difference(coils)
    if len(notfound) > 0:
        logger.info(
            "The requested coil(s) {} are not available ({})".format(notfound, coils)
        )
        return

    # Are there cathcords files in the dst directory? If so, don't overwrite them.
    overlap = set.intersection(set(src_info.keys()), set(dst_info.keys()))
    if len(overlap) > 0:
        logger.info(
            "Found %s already-been-processed recordings. Skipping them.", len(overlap)
        )

    for k in overlap:
        src_info.pop(k, None)

    # - - - - - - - - - - - - -
    # 2. Apply UKF pair-by-pair

    number_of_recordings = len(src_info)

    # TODO Should find a way to pass these params. Catheters may have slightly
    #   different geometry. Also ukf.Q
    ukf = catheter_ukf.UKF()

    with enlighten.get_manager() as manager:

        for j, (index, recording) in enumerate(src_info.items()):
            logger.debug("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
            logger.debug("Processing recording %s", index)

            dist, prox = catheter_utils.cathcoords.read_pair(
                recording[distal], recording[proximal]
            )

            # Extra prep
            obs = numpy.concatenate((dist.coords, prox.coords), axis=1)
            dts = 0.001 * numpy.diff(dist.times, prepend=dist.times[0])

            # Would be nice if there was a way to config this. Initial P especially.
            x, P = ukf.estimate_initial_state(dist.coords[0], prox.coords[0])

            # Setup storage for filter results
            M = dist.coords.shape
            dist_coords_ukf = numpy.zeros(M)
            prox_coords_ukf = numpy.zeros(M)
            number_of_frames = M[0]

            with manager.counter(
                total=number_of_frames,
                desc="Applying UKF {}/{}".format(j, number_of_recordings),
                unit="frames",
                leave=False,
            ) as frame_progress:

                failed = True

                logger.debug("Applying UKF %s/%s ...", j, number_of_recordings)
                for i in range(number_of_frames):
                    try:
                        x, P = ukf.filter(x, P, obs[i], dts[i])
                    except ValueError as err:
                        logger.exception("An exception occured during application of UKF:")
                        logger.warn("Skipping recording number %s", index)
                        break
                    _, dist_coord_ukf, prox_coord_ukf = ukf.tip_and_coils(x)
                    dist_coords_ukf[i, :] = dist_coord_ukf
                    prox_coords_ukf[i, :] = prox_coord_ukf
                    frame_progress.update()
                else:
                    failed = False

                if failed:
                    # just skip this recording
                    continue

                logger.debug("... done applying UKF.")

                # Write the result.
                def _write_result(index, coil, data):
                    filename = catheter_utils.cathcoords.make_filename(index, coil)
                    filename = os.path.join(dst, filename)
                    catheter_utils.cathcoords.write_file(filename, data)

                _write_result(index, distal, dist._replace(coords=dist_coords_ukf))
                _write_result(index, proximal, prox._replace(coords=prox_coords_ukf))
                logger.debug("... finished processing recording %s", index)

            # frame_progress.close()

        logger.info("Done!")


@cathy.command()
@click.argument("path", type=click.Path(exists=True))
@click_log.simple_verbosity_option()
def info(path):
    """Summarize the contents of a projection file."""
    # maybe we want info from another kind of file in the future
    if Path(path).is_dir():
        _proj_dir_info(path)
    else:
        _proj_info(path)


def _proj_dir_info(path):
    discoveries, unknown = catheter_utils.projections.discover_raw(path)
    print("This directory contains:")
    if not discoveries.empty:
        nsri = discoveries["recording"][discoveries["scheme"] == "sri"].nunique()
        if nsri > 0:
            sdf = discoveries[discoveries["scheme"] == "sri"]
            print("  {} 'sri' sampled recordings".format(nsri))
            print("     from coils: {}".format(sdf["coil"].unique()))

        nhadamard = discoveries["recording"][discoveries["scheme"] == "fh"].nunique()
        if nhadamard > 0:
            hdf = discoveries[discoveries["scheme"] == "fh"]
            print("  {} 'fh' sampled recordings".format(nhadamard))
            print("     from coils: {}".format(hdf["coil"].unique()))
            print("   with dithers: {}".format(hdf["dither"].unique()))

    nunknown = len(unknown)
    if nunknown > 0:
        print("  {} unknown .projections files".format(nunknown))

def _proj_info(path):
    # Try to read the file, even if it is corrupt. Lets just see what we can see
    meta, raw, version, corrupt = catheter_utils.projections.read_raw(path, allow_corrupt=True)
    fov = catheter_utils.projections.fov_info(meta, raw)

    if "timestamp" in meta.columns:
        times = meta["timestamp"]
        start, end = (datetime.datetime.fromtimestamp(t//1000) for t in (times.min(), times.max()))
        when = "{} (duration {})".format(start, end-start)
    else:
        when = "unknown date and time"

    print("'legacy_version': {}".format(version))
    print("         corrupt: {}".format(corrupt))
    print("            when: {}".format(when))
    print("             fov: {}".format(fov))
    print("         samples: {}".format(len(raw)))
    print(" readouts/sample: {}".format(len(raw[0])))
    print("  pixels/readout: {}".format(len(raw[0][0])))


@cathy.command()
@click.argument("path", type=click.Path(exists=True))
@click.option("-s", "--scheme", type=click.Choice(["sri", "fh"]), help="Select scheme for visualization.")
@click.option("-c", "--coil", type=Ints(), help="Select coil(s) for visualization.")
@click.option("-a", "--axis", type=Ints(), help="Select axis/axes for visualization.")
@click.option("-d", "--dither", type=Ints(), help="Select dither(s) for visualization")
@click.option("-r", "--recording", type=Ints(), help="Select recording(s) for visualization.")
@click.option("-p", "--pick", type=click.Choice(["rand", "last"]),
              help="How pick a single recording from the final query.")
@click.option("-x", "--xyz", type=click.Path(exists=True))
@click.option("-gt", "--groundtruth", type=click.Path(exists=True), help="Path to ground truth csv.")
@click.option("-en", "--expname", help="str representing experiment name to search ground truth file.")
@click_log.simple_verbosity_option()
def peek(path, pick=None, xyz=None, groundtruth=None, expname=None, **kwargs):
    """Visualize catheter projections.

    If PATH is a directory, select relevant projections from available raw files to visualize.
    If PATH is a .projections file, display the contents of the file.
    """
    if Path(path).is_dir():
        _proj_dir_peek(path, xyz, pick=pick, query=kwargs, gt=groundtruth, exp=expname)
    else:
        _proj_peek(path)


def _proj_dir_peek(path, xyz, *, query, pick, gt=None, exp=None):
    """Peek a .projections directory. Look at available data and select"""
    if pick is None:
        pick = "rand"

    # find all available data in the given directory
    discoveries, unknown = catheter_utils.projections.discover_raw(path)

    # filter data based on user's optional query params
    query = [discoveries[k].isin(v) for k, v in query.items() if v]
    if len(query) == 0:
        selection = discoveries
    elif len(query) == 1:
        selection = discoveries[query[0]]
    else:
        selection = discoveries[functools.reduce(lambda x, y: x & y, query)]

    # narrow selection to a single visualization
    recordings = sorted(selection["recording"].unique())
    if pick == "rand":
        selected_recording = random.choice(recordings)
    else:
        selected_recording = recordings[-1]

    selection = selection[selection["recording"] == selected_recording]

    coils = sorted(selection["coil"].unique())
    axes = sorted(selection["axis"].unique())
    filenames = sorted(selection["filename"].unique())

    if xyz is not None:
        if not Path(xyz).is_dir():
            logger.warning("xyz should be a directory")
            xyz = None

    coord_data = {}
    if xyz is not None:
        try:
            cathcoord_files = catheter_utils.cathcoords.discover_files(xyz)
            for coil in coils:
                coord_data[coil] = catheter_utils.cathcoords.read_file(
                    cathcoord_files[selected_recording][coil]
                )
        except KeyError:
            logger.warning("could not find coordinate file in xyz dir")
            coord_data = {}

    if len(axes) == 3:
        coordinate_system = catheter_utils.localization.XYZCoordinates()
    elif len(axes) == 4:
        coordinate_system = catheter_utils.localization.HadamardCoordinates()
    else:
        logger.warning("Unknown coordinate system")
        coordinate_system = None

    number_of_frames = None
    filename_to_artist = {}
    for filename in filenames:
        meta, raw, _, _ = catheter_utils.projections.read_raw(filename)
        artist = art.ProjectionsAnimator(meta, raw)
        if number_of_frames:
            number_of_frames = min(number_of_frames, len(raw))
        else:
            number_of_frames = len(raw)
        filename_to_artist[filename] = artist

    fig, axs = pyplot.subplots(len(axes), figsize=[9, 6])
    if len(axes) > 1:
        axis_map = {a: axs[i] for i, a in enumerate(axes)}
    else:
        axis_map = {axes[0]: axs}

    def do_nothing(i):
        pass

    def update_axvline(ln1, ln2, coil, proj):
        try:
            coords = coord_data[coil].coords

            def _fn(i):
                pcoord = coordinate_system.world_to_projection(coords[i])
                ln1.set_xdata(pcoord[proj])

                if (gt != None):
                    # search the gt file for experiment name and/or coil index to read ground truth
                    gtcoord = get_gt.read_results(gt, Exp_name=exp, Coil_index=coil) 
                    gtcoord = coordinate_system.world_to_projection(gtcoord)
                    ln2.set_xdata(gtcoord[proj])

            return _fn
        except KeyError:
            return do_nothing

    coil_to_color = {4: 'red', 5: 'blue', 6: 'green', 7: 'purple'}
    cmap =pyplot.cm.get_cmap('tab10')

    plot_keys = []
    for row in selection.itertuples():
        ax = axis_map[row.axis]
        artist = filename_to_artist[row.filename]
        p, = artist.plot(0, row.index, ax=ax, plot_kwargs=dict(color=cmap(row.coil%cmap.N)))

        if coordinate_system is not None and row.coil in coord_data:
            ln1 = ax.axvline(0, color=coil_to_color[row.coil])
            ln2 = ax.axvline(0, color=coil_to_color[row.coil], linestyle='--')
            plot_keys.append((artist, row.index, p, update_axvline(ln1, ln2, row.coil, row.axis)))
        else:
            plot_keys.append((artist, row.index, p, do_nothing))

    def init():
        pass

    def animate(i):
        for artist, index, p, pp in plot_keys:
            artist.animate(i, index, plot=p)
            pp(i)

    a = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=number_of_frames,
        interval=10,
        repeat=True
    )

    pyplot.show()
    pyplot.close(fig)


def _proj_peek(path):
    # Try to read the file, even if it is corrupt. Lets just see what we can see
    meta, raw, version, corrupt = catheter_utils.projections.read_raw(path, allow_corrupt=True)
    if len(raw) == 0 or len(raw[0]) == 0:
        logger.info("no projections found in the file")
        return

    artist = art.ProjectionsAnimator(meta, raw)

    if "timestamp" in meta.columns:
        interval = numpy.mean(numpy.diff(meta["timestamp"]))
    else:
        # if we don't have timestamps assume we are running at approx 24 fps
        # (for 3 readouts). This won't be real if the readouts are stored
        # in multiple files (e.g., the FH case) but what can we do?
        interval = (1000/24) * (len(raw[0])/3)

    fig = pyplot.figure(figsize=[9, 6])
    ax = pyplot.axes(xlim=artist.xlim, ylim=artist.ylim)

    plots = []
    for j in range(len(raw[0])):
        p, = artist.plot(0, j, ax=ax, plot_args=['-*'])
        plots.append(p)

    def init():
        pass

    def animate(i):
        for k, plot in enumerate(plots):
            artist.animate(i, k, plot=plot)

    a = animation.FuncAnimation(fig, animate, init_func=init, frames=len(raw), interval=interval, repeat=True)
    pyplot.show()
    pyplot.close(fig)


@cathy.command()
@click.argument("src_path", type=click.Path(exists=True))
@click.argument("dst_path", type=click.Path())
@click.option("-d", "--distal", "distal_index", type=int, default=5, help="Select distal coil index.")
@click.option("-p", "--proximal", "proximal_index", type=int, default=4, help="Select proximal coil index.")
@click.option("-g", "--geometry", "geometry_index", type=int, default=1, help="Select geometry index.")
@click.option("-z", "--dither", "dither_index", type=int, default=0, help="Select dither index.")
@click_log.simple_verbosity_option()
def localize(src_path, dst_path, distal_index, proximal_index, geometry_index, dither_index):
    toc, unknowns = catheter_utils.projections.discover_raw(src_path)

    dst_path = Path(dst_path)
    if not dst_path.exists():
        dst_path.mkdir()
    if not dst_path.is_dir():
        raise NotADirectoryError()

    geo = catheter_utils.geometry.GEOMETRY[geometry_index]

    # which parameters and where should they come from?
    width = 3.5
    sigma = 0.75

    loc_fns = {
        "peak": _localizer(catheter_utils.localization.peak, None, None),
        "centroid": _localizer(catheter_utils.localization.centroid, None, None),
        "centroid_around_peak": _localizer(catheter_utils.localization.centroid_around_peak, None, dict(window_radius=2*width)),
        "png": _localizer(catheter_utils.localization.png, None, dict(width=width, sigma=sigma)),
        "jpng": _jpng(geo, width=width, sigma=sigma),
        #"wjpng": _wjpng(geo, width=width, sigma=sigma), #unverified
    }

    for recording in sorted(toc.recording.unique()):
        print(f"--- -- -- -- -- -- -- -- -- -- --")
        print(f"[{recording}] pre-processing data")
        data = catheter_utils.projections.FindData(toc, recording, distal_index, proximal_index, dither_index)

        for loc_name, loc_fn in loc_fns.items():
            print(f"[{recording}] processing using {loc_name}")

            distal_coords = []
            proximal_coords = []

            for i in range(len(data)):
                distal_data = data.get_distal_data(i)
                proximal_data = data.get_proximal_data(i)

                distal, proximal = loc_fn(distal_data, proximal_data)
                distal_coords.append(distal)
                proximal_coords.append(proximal)

            distal_coords = numpy.array(distal_coords)
            proximal_coords = numpy.array(proximal_coords)

            for coil, coords in zip(
                    [distal_index, proximal_index],
                    [distal_coords, proximal_coords]
            ):
                cc = catheter_utils.cathcoords.Cathcoords(
                    coords=coords,
                    times=data.timestamp,
                    snr=data.snr,
                    trigs=data.trig,
                    resps=data.resp,
                )

                path = dst_path / Path(loc_name)

                if not path.exists():
                    path.mkdir()

                if not path.is_dir():
                    raise NotADirectoryError()

                fn = catheter_utils.cathcoords.make_filename(recording, coil)
                fn = path / fn

                catheter_utils.cathcoords.write_file(fn, cc)


def _localizer(fn, args, kwargs):
    def _fn(d, p):
        return catheter_utils.localization.localize_catheter(d, p, fn, args, kwargs)
    return _fn


def _jpng(geo, width=3.0, sigma=0.5):
    def _fn(d, p):
        return catheter_utils.localization.jpng(d, p, geo, width=width, sigma=sigma)
    return _fn


def _wjpng(geo, width=3.0, sigma=0.5):
    l = catheter_utils.localization.JointIterativeWeightedCentroid(
        geometry=geo,
        centroid_weighting=(lambda xs: catheter_utils.localization.png_density(xs, width, sigma)),
        err_weighting=catheter_utils.projections.variance,
    )
    return l.localize

@cathy.command()
@click.argument("src_path", type=click.Path(exists=True, file_okay=False))
@click.argument("dicom_file", type=click.Path(exists=True, dir_okay=False))
@click.option("-d", "--distal", "distal_index", type=int, default=5, help="Select distal coil index.")
@click.option("-p", "--proximal", "proximal_index", type=int, default=4, help="Select proximal coil index.")
@click.option("-r", "--recording", "recording_index", type=int, help="Select recording index.")
@click.option("-g", "--geometry", "geometry_index", type=int, default=1, help="Select geometry index.")
@click.option("-gt", "--groundtruth", type=click.Path(exists=True), default=None)
@click.option("-en", "--expname", default=None)
@click_log.simple_verbosity_option()
def scatter(src_path, dicom_file, distal_index, proximal_index, recording_index, geometry_index, groundtruth, expname):

    cathcoord_files = catheter_utils.cathcoords.discover_files(src_path)
    if recording_index is None:
        recording_index = random.choice(list(cathcoord_files.keys()))

    distal_file = cathcoord_files[recording_index][distal_index]
    proximal_file = cathcoord_files[recording_index][proximal_index]

    distal, proximal = catheter_utils.cathcoords.read_pair(distal_file, proximal_file)
    geo = catheter_utils.geometry.GEOMETRY[geometry_index]
    fit = geo.fit_from_coils_mse(distal.coords, proximal.coords)

    if groundtruth != None: 
        distal_gtcoord = get_gt.read_results(groundtruth, Exp_name=expname, Coil_index=distal_index)
        proximal_gtcoord = get_gt.read_results(groundtruth, Exp_name=expname, Coil_index=proximal_index)
        fit_gt = geo.fit_from_coils_mse(distal_gtcoord, proximal_gtcoord)

    d = dicom.read_file(dicom_file)
    art = dicom_art.DicomPlotter(d)

    pyplot.figure()
    art.imshow()
    art.plot(fit.tip, plot_args=[".g"])
    art.plot(fit.distal, plot_args=[".b"])
    art.plot(fit.proximal, plot_args=[".r"])

    if groundtruth != None:
        art.plot(fit_gt.tip, plot_args=["*g"])
        art.plot(fit_gt.distal, plot_args=["*b"])
        art.plot(fit_gt.proximal, plot_args=["*r"])
    
    pyplot.show()


@cathy.command() 
@click.argument("args", nargs=-1)
@click.option("-d", "--dest", type=click.Path(exists=True, file_okay=False), default="./")
@click.option("-en", "--expname", default=None)
@click.option("-m", "--meta", default="")
@click.option("-c", "--coil", type=int, default=None)
@click_log.simple_verbosity_option()
def gt_tool(args, dest, expname, meta, coil):
    get_gt.final_coords(args, dest_path=dest, Exp_name=expname, Meta_data=meta, Coil_index=coil)


@cathy.command() 
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click_log.simple_verbosity_option()
def view_dicom(path):
    get_gt.quick_view(path)

@cathy.group()
def coil_metrics():
    """Coil Metrics subcommands"""
    pass
@coil_metrics.command()
@click.argument("tracking_sequence_folder", type=click.Path(exists=True))
@click.argument("localization_folder", type=click.Path(exists=True))
@click.argument("groundtruth_folder", type=click.Path(exists=True))
@click.argument("algorithms")
@click.option("--dest", type=click.Path(exists=True, file_okay=False), help="output directory path for trackingerr.csv file. If not specified default is localization algorithm results folder", default=None)
@click.option("--distal_index", "-d", default=5, help="Distal coil index.")
@click.option("--proximal_index", "-p", default=4, help="Proximal coil index.")
@click.option("-z", "--dither", "dither_index", type=int, default=0, help="Select dither index.")
@click.option("-en", "--expname", default=None, help="experiment name in GroundTruthCoords.csv")

@click_log.simple_verbosity_option()
def calc(tracking_sequence_folder, localization_folder, groundtruth_folder, algorithms, dest, distal_index, proximal_index, dither_index, expname):
    '''
    calculate or plot the tracking error of localization algorithms based on bias and 95% ChebyShev (compare stats of groundtruth to localization results)

    \b
    tracking_sequence_folder -> path to tracking sequence folder containing .projection files
    localization_folder -> path to localization algorithm results folder
    groundtruth_folder -> path to folder with GroundTruthCoords.csv ... add "/" to end of path (PATH/GroundTruthFolder/)
    algorithms -> localization algorithms EX: peak,centroid_around_peak,png .Note: Do not use any spaces between commas
    '''

    catheter_utils.metrics.trackerr(tracking_sequence_folder, localization_folder, groundtruth_folder, algorithms, dest, distal_index, proximal_index, dither_index, expname)

@coil_metrics.command()
@click.argument("trackerr_folder")
@click.option("-en", "--expname", default=None, help="experiment name in GroundTruthCoords.csv")
@click.option("-x", "--xaxis", default='algorithm', type=click.Choice(['algorithm', 'FOV', 'dither', 'trackseq']), help="choose x-axis column for plotting. Options: algorithm, FOV, dither, or trackseq")
@click_log.simple_verbosity_option()
def plotter(trackerr_folder, expname, xaxis):
    '''
    Barplot of tracking errors

    trackerr_folder -> path to folder containing *trackingerr.csv file
    '''

    catheter_utils.metrics.barplot(trackerr_folder, expname, xaxis)

@cathy.command()
@click.argument("tracking_sequence_folder_path", type=click.Path(exists=True))
@click.argument("dest", type=click.Path(exists=True))
@click.argument("recording_indexes", type=Ints())
@click.option("-d", "--distal", "distal_index", type=int, default=5, help="Select distal coil index.")
@click.option("-p", "--proximal", "proximal_index", type=int, default=4, help="Select proximal coil index.")
@click.option("-z", "--dither", "dither_index", type=int, default=None, help="Select dither index.")

@click_log.simple_verbosity_option()
#TODO select and manage multiple dither indexes/folders
def fft_projections(tracking_sequence_folder_path, dest, recording_indexes, distal_index, proximal_index, dither_index):
    '''Get absolute fft signal data from projection files for specified dither and distal & proximal coils and export as text file into fft subfolder

    \b
    Textfile name template is "fft_signals_coilX_recX_fovXXX_readoutsXXX.txt" (coil index, recording index, field of view (mm), # of readouts)

    Example of text file format
    PRJTXT (file format name)
    1 (version)
    SRI_April-2020-07-09T12_20_56.765 (tracking sequence folder name)
    x y z x_timestamp average_timestamp readout_size (column headings: axes... length of each readout)
    '''
    dithered = 0
    subfolder = "fft"
    # create subfolder and specify by dither index (if applicable)
    if isinstance(dither_index, int):
        dithered = dither_index
        subfolder = "/".join((subfolder, "dith" + str(dithered)))
    #add subfolder to destination path
    dest = os.path.join(dest, subfolder)
    if not os.path.exists(dest):
        os.makedirs(dest)
    if not os.path.isdir(dest):
        raise Exception("Error with saving to output directory. Check output directory path")

    # find projection files
    discoveries, _ = catheter_utils.projections.discover_raw(tracking_sequence_folder_path)
    # check that specified coils, recording & dither exist in testdata set
    assert set({distal_index, proximal_index}).issubset(set(discoveries["coil"].values)) and dithered in discoveries['dither'].values and set(recording_indexes).issubset(set(discoveries['recording'].values)), \
        f"specified distal/proximal coils {distal_index, proximal_index}. coils found are {set(discoveries['coil'].unique())}. " \
        f"Recording specified {recording_indexes}. Recording(s) found {set(discoveries['recording'].unique())}. Dither specified {dithered}. Dither(s) found {set(discoveries['dither'].unique())}"

    # extract projections
    # FindData -> organize projection data & perform absolute fourier transform/shift
    for record in recording_indexes:
        data = catheter_utils.projections.FindData(discoveries, record, distal_index, proximal_index, dithered)
        #extract proximal and distal data
        pmeta, _, _ = data._data[(proximal_index, 0)]
        dmeta, _, _ = data._data[(distal_index, 0)]
        _export_data(dest, tracking_sequence_folder_path.split('/')[-1], distal_index,  "distal", record, pmeta.fov[0], data)
        _export_data(dest, tracking_sequence_folder_path.split('/')[-1], proximal_index, "proximal", record, dmeta.fov[0], data)

def _export_data(dir, seqname, coil, coil_name, rec, fov, prj):
    '''
    export absolute magnitude fft projection signal data to text file
    '''
    fs = numpy.array([])
    columns = ['x', 'y', 'z',  'x-timestamp', 'average-timestamp', 'readout_size']
    # format and version
    format = "PRJTXT"
    version = 1
    # extract all readout details to an array
    for readout in range(len(prj)):
        data = prj._get_coil_data(coil, readout, getattr(prj, f"_{coil_name}_cache"))
        coords = numpy.array(data)[:, 0, :].T
        coords_len = numpy.tile(len(coords), (len(coords), 1))
        timestamp_first= numpy.tile(int(prj.get_timestamp(0,readout)), (len(coords), 1))
        timestamp_ave = numpy.tile(int(prj._timestamp[readout]), (len(coords), 1))
        if readout == 0:
            fs = numpy.column_stack((coords, timestamp_first, timestamp_ave ,coords_len))
        else:
            fs = numpy.vstack((fs, numpy.column_stack((coords, timestamp_first, timestamp_ave ,coords_len))))

    # create data table and export as textfile
    df = pd.DataFrame(fs)
    if fs.shape[1] > 6: # for hadamard add s axis to header [x, y, z, s, timestamp, timestamp average, readout_size]
        columns.insert(3, 's')
    filefft = os.path.join(dir, f"fft_signals_coil{coil}_rec{rec}_fov{int(fov)}_readouts{len(prj)}.txt")

    with open(filefft, "w") as f:
        f.write(f"{format}\n")
        f.write(f"{version}\n")
        f.write(f"{seqname}\n")
        f.write(f"{' '.join(columns)}\n")
        df.to_csv(f, header=False, sep=" ", index=False)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cathy build-resp

# @cathy.command()
# def build_resp():
#     pass


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# cathy apply-resp

# @cathy.command()
# def apply_resp():
#     pass
