"""Command Line Interface for cathy."""
import collections
import datetime
import functools
import logging
import os
import random
import re
from functools import reduce
from pathlib import Path
from pprint import pprint

import click
import click_log
import enlighten
import numpy
import pandas
from matplotlib import animation, pyplot

import catheter_ukf
import catheter_utils.cathcoords
import catheter_utils.geometry
import catheter_utils.localization
import catheter_utils.projections
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

    if "fov" in meta.columns:
        pix = len(raw[0][0])
        if meta["fov"].nunique() == 1:
            fov = meta["fov"][0]
            fov = "{} mm (resolution {:.2f} mm)".format(fov, fov/pix)
        else:
            min_fov, max_fov = meta["fov"].min(), meta["fov"].max()
            fov = "{} - {} (resolution {:.2f} - {:.2f} mm)".format(min_fov, max_fov, min_fov/pix, max_fov/pix)
    else:
        fov = "unknown fov"

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


# def _select_cath_data(toc, coils, )


@cathy.command()
@click.argument("path", type=click.Path(exists=True))
@click.option("-s", "--scheme", type=click.Choice(["sri", "fh"]), help="Select scheme for visualization.")
@click.option("-c", "--coil", type=Ints(), help="Select coil(s) for visualization.")
@click.option("-a", "--axis", type=Ints(), help="Select axis/axes for visualization.")
@click.option("-d", "--dither", type=Ints(), help="Select dither(s) for visualization")
@click.option("-r", "--recording", type=Ints(), help="Select recording(s) for visualization.")
@click.option("-p", "--pick", type=click.Choice(["rand", "last"]),
              help="How pick a single recording from the final query.")
@click_log.simple_verbosity_option()
def peek(path, pick=None, **kwargs):
    """Visualize catheter projections.

    If PATH is a directory, select relevant projections from available raw files to visualize.
    If PATH is a .projections file, display the contents of the file.
    """
    if Path(path).is_dir():
        _proj_dir_peek(path, pick=pick, query=kwargs)
    else:
        _proj_peek(path)


def _proj_dir_peek(path, *, query, pick):
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

    plot_keys = []
    for row in selection.itertuples():
        ax = axis_map[row.axis]
        artist = filename_to_artist[row.filename]
        p, = artist.plot(0, row.index, ax=ax)
        plot_keys.append((artist, row.index, p))

    def init():
        pass

    def animate(i):
        for artist, index, p in plot_keys:
            artist.animate(i, index, plot=p)

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
@click.argument("path", type=click.Path(exists=True))
@click_log.simple_verbosity_option()
def localize(path):
    toc, unknowns = catheter_utils.projections.discover_raw(path)

    # filter using a query
    distal_index = 4
    proximal_index = 5
    coil_indices = {distal_index, proximal_index}

    # separable localization algorithms
    loc_fns = {
        "peak": catheter_utils.localization.peak,
        "centroid": catheter_utils.localization.centroid,
        "centroid_around_peak": catheter_utils.localization.centroid_around_peak,
        "png": catheter_utils.localization.png,
    }

    for recording in sorted(toc.recording.unique()):
        if recording != 0:
            continue

        data = _read_recording(toc, recording, coil_indices)

        axes = set()
        loc_data = collections.defaultdict(list)
        for (coil, axis, dither), (meta, projections) in data.items():
            axes.add(axis)
            for (row, signal) in zip(meta.itertuples(), projections):
                xs = numpy.linspace(-row.fov / 2, row.fov / 2, len(signal))
                for loc, fn in loc_fns.items():
                    x = fn(signal, xs)
                    loc_data[(coil, axis, dither, loc)].append(x)

        for k, v in loc_data.items():
            loc_data[k] = numpy.array(v)

        for loc in loc_fns.keys():
            for coil in coil_indices:
                a = []
                for axis in range(3):
                    meta, _ = data[(coil, axis, 0)]

                    timestamp, resp, trig = _get_column_or_zeros(
                        meta,
                        ["timestamp", "resp", "trig"]
                    )

                    a.append(pandas.DataFrame({
                        "timestamp": timestamp,
                        f"c{axis}": loc_data[(coil, axis, 0, loc)],
                        f"resp{axis}": resp,
                        f"trig{axis}": trig,
                    }))

                df = reduce(lambda u, v: pandas.merge(u, v, on="timestamp"), a)
                times = df["timestamp"].values
                resp = numpy.mean(df[[f"resp{a}" for a in axes]].values, axis=1)
                trig = numpy.mean(df[[f"trig{a}" for a in axes]].values, axis=1)

                cc = catheter_utils.cathcoords.Cathcoords(
                    coords=df[["c0", "c1", "c2"]].values,
                    times=times,
                    snr=numpy.zeros(len(times)),
                    trigs=trig,
                    resps=resp
                )

                path = Path(loc)

                if not path.exists():
                    path.mkdir()

                if not path.is_dir():
                    raise NotADirectoryError()

                fn = catheter_utils.cathcoords.make_filename(recording, coil)
                fn = path / fn

                catheter_utils.cathcoords.write_file(fn, cc)


def _get_column_or_zeros(df, names):
    res = []
    n = len(df.index)
    for name in names:
        if name in df:
            res.append(df[name].values)
        else:
            res.append(numpy.zeros(n))
    return res


def _read_raw_and_reconstruct(fn):
    meta, projections, _, _ = catheter_utils.projections.read_raw(fn, allow_corrupt=False)
    for i, signal in enumerate(projections):
        projections[i] = catheter_utils.projections.reconstruct(projections[i])
    return meta, projections


def _read_recording(toc, recording, coil_indices):
    """Read recording data related to the given coils."""

    selection = toc[(toc.recording == recording) & toc.coil.isin(coil_indices)]

    filename_to_data = {
        fn: _read_raw_and_reconstruct(fn)
        for fn in sorted(selection.filename.unique())
    }

    data = {}
    for row in selection.itertuples():
        meta, projections = filename_to_data[row.filename]
        sliced_projections = [p[row.index] for p in projections]
        data[(row.coil, row.axis, row.dither)] = meta, sliced_projections

    return data

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
