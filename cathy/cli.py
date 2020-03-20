"""Command Line Interface for cathy."""

import catheter_ukf
import catheter_utils.cathcoords
import catheter_utils.projections
import click
import click_log
import datetime
import enlighten
import logging
import numpy
import os

from matplotlib import animation
from matplotlib import pyplot


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
    # different geometry. Also ukf.Q
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
    _proj_info(path)


def _proj_info(path):
    # Try to read the file, even if it is corrupt. Lets just see what we can see
    meta, raw, version, corrupt = catheter_utils.projections.read_raw(path, allow_corrupt=True)

    if "fov" in meta.columns:
        pix = len(raw[0][0])
        if meta["fov"].nunique() == 1:
            fov = meta["fov"][0]
            fov = "{} mm (resolution {} mm)".format(fov, fov/pix)
        else:
            min_fov, max_fov = meta["fov"].min(), meta["fov"].max()
            fov = "{} - {} (resolution {} - {} mm)".format(min_fov, max_fov, min_fov/pix, max_fov/pix)
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


@cathy.command()
@click.argument("path", type=click.Path(exists=True))
@click_log.simple_verbosity_option()
def peek(path):
    """Display the contents of a projection file."""
    # maybe we want to peek into another kind of file in the future
    _proj_peek(path)


def _proj_peek(path):
    # Try to read the file, even if it is corrupt. Lets just see what we can see
    meta, raw, version, corrupt = catheter_utils.projections.read_raw(path, allow_corrupt=True)
    if len(raw) == 0 or len(raw[0]) == 0:
        logger.info("no projections found in the file")
        return

    def pick(i):
        fov = meta["fov"][i]
        fs = catheter_utils.projections.reconstruct(raw[i])
        xs = numpy.linspace(-fov / 2, fov / 2, len(fs[0]))
        return fs, xs

    f0, x0 = pick(0)
    if "timestamp" in meta.columns:
        interval = numpy.mean(numpy.diff(meta["timestamp"]))
    else:
        # if we don't have timestamps assume we are running at approx 24 fps
        # (for 3 readouts). This won't be real if the readouts are stored
        # in multiple files (e.g., the FH case) but what can we do?
        interval = (1000/24) * (len(f0)/3)

    fig = pyplot.figure(figsize=[15, 10])
    ax = pyplot.axes(xlim=(x0[0], x0[-1]), ylim=(0, numpy.max(f0)))

    plots = []
    for j in range(len(f0)):
        p, = ax.plot(x0, f0[j], '-*')
        plots.append(p)

    def init():
        pass

    def animate(i):
        fs, xs = pick(i)
        for k, plot in enumerate(plots):
            plot.set_data(xs, fs[k])

    a = animation.FuncAnimation(fig, animate, init_func=init, frames=len(raw), interval=interval, repeat=True)
    pyplot.show()
    pyplot.close(fig)

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
