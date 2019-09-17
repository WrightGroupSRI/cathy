"""Command line tools for manipulating catheter files."""

import catheter_ukf
import catheter_utils
import catheter_uniform_rmm
import dicom_utils

import matplotlib.pyplot as plt
import matplotlib.animation
import numpy
import pydicom as dicom
import traceback
import seaborn as sns
import scipy.stats

from collections import defaultdict


# Hard coded dumb stuff
root = "/media/snferguson/Data Volume/Pigs/MR-EP-TC-2017-09-23-GIP57/Visualization"
vmap_dirname = root + "/voltageMap/SRI_Catheter_Tracking-2017-09-23T15:36:13.785"
abln_dirname = root + "/ablation1/SRI_Catheter_Tracking-2017-09-23T14:22:59.237"
sax_cine_dirname = root + "/Volumes/E3891S010"
lax_cine_dirname = root + "/Volumes/E3891S011"
output_dirname = "/home/snferguson/Desktop/movies"
cache_dir = output_dirname + "/cache"

dist_coil = 4  # tip adjacent
prox_coil = 5
trig_window = 375, 425
ukf = catheter_ukf.UKF()


def get_nearest_key(d, k):
    return d[min(d.keys(), key=lambda u: abs(u - k))]


def calculate_filtered_coords(ukf, dist_coords, prox_coords, dts):
    # init filter state
    x, P = ukf.estimate_initial_state(dist_coords[0], prox_coords[0])

    # setup storage for filter results
    M = dist_coords.shape
    dist_coords_ukf = numpy.zeros(M)
    prox_coords_ukf = numpy.zeros(M)
    tip_coords_ukf = numpy.zeros(M)

    # run ukf over all the data
    for i, (dist_coord, prox_coord, dt) in enumerate(
        zip(dist_coords, prox_coords, dts)
    ):
        obs = numpy.concatenate((dist_coord, prox_coord))
        x, P = ukf.filter(x, P, obs, dt)

        tip_coord_ukf, dist_coord_ukf, prox_coord_ukf = ukf.tip_and_coils(x)
        dist_coords_ukf[i, :] = dist_coord_ukf
        prox_coords_ukf[i, :] = prox_coord_ukf
        tip_coords_ukf[i, :] = tip_coord_ukf

    return tip_coords_ukf, dist_coords_ukf, prox_coords_ukf


CCOLOR = "#fc00ef"
UCOLOR = "#00fcfc"
MCOLOR = "#fc9f00"


class Plotter(object):
    def __init__(
        self,
        tip_coords,
        dist_coords,
        prox_coords,
        tip_coords_ukf,
        dist_coord_ukf,
        prox_coord_ukf,
        triggers,
        cine_data,
        label="unfiltered",
        label_ukf="UKF",
    ):
        self.tip_coords = tip_coords
        self.dist_coords = dist_coords
        self.prox_coords = prox_coords
        self.tip_coords_ukf = tip_coords_ukf
        self.dist_coords_ukf = dist_coords_ukf
        self.prox_coords_ukf = prox_coords_ukf
        self.trigs = triggers

        # pick image slice nearest to the points
        centroid = numpy.mean(self.tip_coords_ukf, axis=0)
        exemplar = list(list(cine_data.values())[0].values())[0]
        self.slice_data = get_nearest_key(
            cine_data, dicom_slice_direction(exemplar, centroid)
        )
        self.label = label
        self.label_ukf = label_ukf

    def _draw_cath(self, tip_coord, dist_coord, prox_coord, color, label):
        plt.plot(
            [prox_coord[0], tip_coord[0]],
            [-prox_coord[1], -tip_coord[1]],
            "-",
            color=color,
            label=label,
        )
        plt.plot(tip_coord[0], -tip_coord[1], "x", color=color)
        plt.plot(dist_coord[0], -dist_coord[1], ".", color=color)
        plt.plot(prox_coord[0], -prox_coord[1], ".", color=color)

    def _draw_instant(self, n):

        tip_coord, dist_coord, prox_coord = (
            self.tip_coords[n],
            self.dist_coords[n],
            self.prox_coords[n],
        )
        tip_coord_ukf, dist_coord_ukf, prox_coord_ukf = (
            self.tip_coords_ukf[n],
            self.dist_coords_ukf[n],
            self.prox_coords_ukf[n],
        )

        image = get_nearest_key(self.slice_data, self.trigs[n])
        image_transform = numpy.linalg.inv(dicom_transform(image))

        tip_coord = image_transform @ numpy.array([*tip_coord, 1])
        dist_coord = image_transform @ numpy.array([*dist_coord, 1])
        prox_coord = image_transform @ numpy.array([*prox_coord, 1])
        tip_coord_ukf = image_transform @ numpy.array([*tip_coord_ukf, 1])
        dist_coord_ukf = image_transform @ numpy.array([*dist_coord_ukf, 1])
        prox_coord_ukf = image_transform @ numpy.array([*prox_coord_ukf, 1])

        plt.imshow(image.pixel_array, cmap="gray")
        self._draw_cath(
            tip_coord_ukf, dist_coord_ukf, prox_coord_ukf, UCOLOR, label=self.label_ukf
        )
        self._draw_cath(tip_coord, dist_coord, prox_coord, CCOLOR, label=self.label)

    def _update(self, n):
        plt.cla()
        # plt.axis((140, 370, 182, 365))
        self._draw_instant(n)
        plt.legend()

    def render_to_file(self, name, *args, **kwargs):
        fig, _ = plt.subplots()
        mw = matplotlib.animation.AVConvFileWriter(fps=kwargs.get("fps", 23))
        with mw.saving(fig, name, dpi=kwargs.get("dpi", 200)):
            for n in range(len(self.tip_coords)):
                self._update(n)
                mw.grab_frame()


class Plotter4x4(object):
    def __init__(
        self,
        cath_data,
        cath_data_ukf,
        cath_data_resp,
        cath_data_ukf_resp,
        triggers,
        sax_cine_data,
        lax_cine_data,
    ):
        self.cath_plotter_sax = Plotter(
            *cath_data, *cath_data_ukf, triggers, sax_cine_data
        )
        self.cath_plotter_lax = Plotter(
            *cath_data, *cath_data_ukf, triggers, lax_cine_data
        )
        self.cath_plotter_sax_resp = Plotter(
            *cath_data_resp, *cath_data_ukf_resp, triggers, sax_cine_data
        )
        self.cath_plotter_lax_resp = Plotter(
            *cath_data_resp, *cath_data_ukf_resp, triggers, lax_cine_data
        )
        self.N = len(triggers)

    def render_to_file(self, name, *args, **kwargs):
        fig, _ = plt.subplots()

        # fig, ax = plt.subplots(2, 2, sharex='col', sharey='row')

        mw = matplotlib.animation.ImageMagickWriter(fps=kwargs.get("fps", 23))
        with mw.saving(fig, name, dpi=kwargs.get("dpi", 200)):
            for n in range(self.N):
                plt.subplot(2, 2, 1)
                self.cath_plotter_sax._update(n)

                plt.subplot(2, 2, 2)
                self.cath_plotter_lax._update(n)

                plt.subplot(2, 2, 3)
                self.cath_plotter_sax_resp._update(n)

                plt.subplot(2, 2, 4)
                self.cath_plotter_lax_resp._update(n)

                mw.grab_frame()


vmap_files = discover_cathcoords_files(vmap_dirname)
abln_files = discover_cathcoords_files(abln_dirname)
join_files = vmap_files.copy()
join_files[9999] = abln_files[0]


sax_cine_data = read_cine_dir(sax_cine_dirname)
lax_cine_data = read_cine_dir(lax_cine_dirname)


bad_tracking = set([14, 22, 23, 25, 43, 44])
bad_triggering = set([11, 47, 70, 95, 124, 146])
moved = set([26])
bad_data = bad_tracking.union(moved, bad_triggering)


def data_for_set(key, join_files):
    name = (cache_dir + "/processed-{:04d}.pickle").format(key)
    try:
        with open(name, "rb") as f:
            tip_coords, dist_coords, prox_coords, tip_coords_ukf, dist_coords_ukf, prox_coords_ukf, dts, trigs, resps = pickle.load(
                f
            )
            print("Read cached data for set {}".format(key))
            return (
                tip_coords,
                dist_coords,
                prox_coords,
                tip_coords_ukf,
                dist_coords_ukf,
                prox_coords_ukf,
                dts,
                trigs,
                resps,
            )
    except FileNotFoundError:
        pass
    except:
        print("An exception occured while trying to load!")
        traceback.print_exc()

    print(
        "Reading files:\n  {}\n  {}".format(
            join_files[key][dist_coil], join_files[key][prox_coil]
        )
    )
    dist_coords, prox_coords, dts, trigs, resps = read_cathcoords_data(
        join_files[key][dist_coil], join_files[key][prox_coil]
    )

    print("Estimating tip the old way")
    M = dist_coords.shape
    tip_coords = numpy.zeros(M)
    for i in range(M[0]):
        delta = dist_coords[i] - prox_coords[i]
        delta = delta / numpy.linalg.norm(delta, ord=2)
        tip_coords[i] = (
            dist_coords[i]
            + (ukf._ukf.statespace.tip_offset - ukf._ukf.statespace.coil_offset) * delta
        )

    print("Applying UKF")
    tip_coords_ukf, dist_coords_ukf, prox_coords_ukf = calculate_filtered_coords(
        ukf, dist_coords, prox_coords, dts
    )

    x = (
        tip_coords,
        dist_coords,
        prox_coords,
        tip_coords_ukf,
        dist_coords_ukf,
        prox_coords_ukf,
        dts,
        trigs,
        resps,
    )

    with open(name, "wb") as f:
        pickle.dump(x, f)

    return (
        tip_coords,
        dist_coords,
        prox_coords,
        tip_coords_ukf,
        dist_coords_ukf,
        prox_coords_ukf,
        dts,
        trigs,
        resps,
    )


# def thing():
#     print("building motion model")


#     coords = None
#     resps = None
#     trigs = None
#     sets = None

#     for key in vmap_files.keys():
#         try:
#             if key in bad_data:
#                 continue
#             _, _, _, c, _, _, _, t, r = data_for_set(key, vmap_files)
#             s = numpy.repeat(key, len(r))
#             coords = numpy.concatenate((coords, c)) if coords is not None else c
#             resps = numpy.concatenate((resps, r)) if resps is not None else r
#             trigs = numpy.concatenate((trigs, t)) if trigs is not None else t
#             sets = numpy.concatenate((sets, s)) if sets is not None else s
#         except:
#             continue

#     # TODO apply gating.
#     print("Gating points")
#     inds = numpy.logical_and(trigs >= trig_window[0], trigs <= trig_window[1])
#     coords = coords[inds]
#     resps = resps[inds]
#     sets = sets[inds]

#     center = numpy.mean(coords, axis=0)
#     print("center = {}".format(center))

#     print("building model")
#     oparams, _, _, _ = pbm.make_model(coords, resps, sets, center=center)
#     print("model params = {}".format(oparams))

#     print("model built")
#     return oparams

# model_params = thing()
# print("Model params:")
# print(model_params)
# quit()


pbm.model.extra_params["center"] = numpy.array(
    [-33.31244365, -98.25879082, -55.93848791]
)
model_params = numpy.array(
    [
        9.99297567e-01,
        1.90996988e-04,
        2.13303515e-02,
        3.08115741e-02,
        2.93206727e00,
        1.44883265e00,
        2.07494068e00,
    ]
)
min_phase = -28876
max_phase = 10586


# slice_index = 10
# slice_data = sax_cine_data[list(sax_cine_data.keys())[slice_index]]
# for key in sorted(slice_data.keys()):
#     name = "/000000000-sax-{}.png".format(int(key))
#     plt.imshow(slice_data[key].pixel_array, cmap="gray")
#     plt.savefig(output_dirname+name)
# quit()

# frame_index = 20
# if True:


# std_devs = None
# std_devs_ukf = None
# std_devs_resp = None
# std_devs_ukf_resps = None


# errs = None
# errs_ukf = None
# errs_resp = None
# errs_ukf_resps = None
# errs_domain = None

# # frame_index = 9999
# # if True:
# for frame_index in join_files.keys():
#     if frame_index in bad_data:
#         continue
#     if frame_index == 9999:
#         continue

#     # print("Processing frame {}.".format(frame_index))

#     try:
#         tip_coords, dist_coords, prox_coords, tip_coords_ukf, dist_coords_ukf, prox_coords_ukf, dts, trigs, resps = data_for_set(frame_index, join_files)
#         resps = pbm.normalize_phase(resps, min_phase=min_phase, max_phase=max_phase)

#         tip_coords_resp = pbm.transform_xyz(tip_coords.copy(), resps, model_params)
#         dist_coords_resp = pbm.transform_xyz(dist_coords.copy(), resps, model_params)
#         prox_coords_resp = pbm.transform_xyz(prox_coords.copy(), resps, model_params)

#         tip_coords_ukf_resp = pbm.transform_xyz(tip_coords_ukf.copy(), resps, model_params)
#         dist_coords_ukf_resp = pbm.transform_xyz(dist_coords_ukf.copy(), resps, model_params)
#         prox_coords_ukf_resp = pbm.transform_xyz(prox_coords_ukf.copy(), resps, model_params)

#         inds = numpy.logical_and(trigs >= trig_window[0], trigs <= trig_window[1])
#         trigs = trigs[inds]
#         resps = resps[inds]

#         tip_coords = tip_coords[inds]
#         dist_coords = dist_coords[inds]
#         prox_coords = prox_coords[inds]

#         tip_coords_resp = tip_coords_resp[inds]
#         dist_coords_resp = dist_coords_resp[inds]
#         prox_coords_resp = prox_coords_resp[inds]

#         tip_coords_ukf = tip_coords_ukf[inds]
#         dist_coords_ukf = dist_coords_ukf[inds]
#         prox_coords_ukf = prox_coords_ukf[inds]

#         tip_coords_ukf_resp = tip_coords_ukf_resp[inds]
#         dist_coords_ukf_resp = dist_coords_ukf_resp[inds]
#         prox_coords_ukf_resp = prox_coords_ukf_resp[inds]


#         mean_err = numpy.mean(tip_coords, axis=0)
#         err = numpy.linalg.norm(tip_coords - mean_err, ord=2, axis=1)**2
#         errs = numpy.concatenate((errs, err)) if errs is not None else err

#         mean_err_ukf = numpy.mean(tip_coords_ukf, axis=0)
#         err_ukf = numpy.linalg.norm(tip_coords_ukf - mean_err_ukf, ord=2, axis=1)**2
#         errs_ukf = numpy.concatenate((errs_ukf, err_ukf)) if errs_ukf is not None else err_ukf

#         mean_err_resp = numpy.mean(tip_coords_resp, axis=0)
#         err_resp = numpy.linalg.norm(tip_coords_resp - mean_err_resp, ord=2, axis=1)**2
#         errs_resp = numpy.concatenate((errs_resp, err_resp)) if errs_resp is not None else err_resp

#         mean_err_ukf_resp = numpy.mean(tip_coords_ukf_resp, axis=0)
#         err_ukf_resp = numpy.linalg.norm(tip_coords_ukf_resp - mean_err_ukf_resp, ord=2, axis=1)**2
#         errs_ukf_resps = numpy.concatenate((errs_ukf_resps, err_ukf_resp)) if errs_ukf_resps is not None else err_ukf_resp

#         errs_domain = numpy.concatenate((errs_domain, resps)) if errs_domain is not None else resps


#         std_dev = numpy.array([numpy.mean(err)**0.5])
#         std_devs = numpy.concatenate((std_devs, std_dev)) if std_devs is not None else std_dev

#         std_dev_ukf = numpy.array([numpy.mean(err_ukf)**0.5])
#         std_devs_ukf = numpy.concatenate((std_devs_ukf, std_dev_ukf)) if std_devs_ukf is not None else std_dev_ukf

#         std_dev_resp = numpy.array([numpy.mean(err_resp)**0.5])
#         std_devs_resp = numpy.concatenate((std_devs_resp, std_dev_resp)) if std_devs_resp is not None else std_dev_resp

#         std_dev_ukf_resp = numpy.array([numpy.mean(err_ukf_resp)**0.5])
#         std_devs_ukf_resps = numpy.concatenate((std_devs_ukf_resps, std_dev_ukf_resp)) if std_devs_ukf_resps is not None else std_dev_ukf_resp


#     except KeyboardInterrupt:
#         raise

#     except:
#         print("Something went wrong ... moving to the next frame")
#         traceback.print_exc()
#         continue


# fig, ax = plt.subplots(2, 2, sharey='all')
# # plt.subplot(2, 2, 1)
# ax[0, 0].plot(errs_domain, errs, '.')
# ax[0, 0].title.set_text("No Correction")
# # plt.subplot(2, 2, 2)
# ax[0, 1].plot(errs_domain, errs_ukf, '.')
# ax[0, 1].title.set_text("UKF")
# # plt.subplot(2, 2, 3)
# ax[1, 0].plot(errs_domain, errs_resp, '.')
# ax[1, 0].title.set_text("Resp")
# # plt.subplot(2, 2, 4)
# ax[1, 1].plot(errs_domain, errs_ukf_resps, '.')
# ax[1, 1].title.set_text("UKF + Resp")
# plt.show()


# df = pandas.DataFrame(data={"No Correction": std_devs, "UKF": std_devs_ukf, "Resp": std_devs_resp, "UKF + Resp": std_devs_ukf_resps}, index=range(len(std_devs)))

# print("- - - - - - - - - -")

# print(scipy.stats.ttest_rel(std_devs, std_devs_ukf))

# print(scipy.stats.ttest_rel(std_devs, std_devs_resp))

# print(scipy.stats.ttest_rel(std_devs_ukf, std_devs_ukf_resps))

# print(scipy.stats.ttest_rel(std_devs_resp, std_devs_ukf_resps))

# # qs = [0.95, 0.9, 0.75, 0.5, 0.25]
# # print(df["No Correction"].std())
# # print(df["UKF"].std())
# # print(df["Resp"].std())
# # print(df["UKF + Resp"].std())

# fix, ax = plt.subplots()
# # df.boxplot(ax=ax)
# # sns.violinplot(data=df)
# # plt.ylabel = "Average Error (mm)"
# # plt.xlabel = "Correction Method"
# ax.set_title("Standard deviation (mm) vs Correction Method")
# plt.show()


for frame_index in join_files.keys():
    if frame_index in bad_data:
        continue

    print("Processing frame {}.".format(frame_index))

    try:
        tip_coords, dist_coords, prox_coords, tip_coords_ukf, dist_coords_ukf, prox_coords_ukf, dts, trigs, resps = data_for_set(
            frame_index, join_files
        )

        print("Applying respiratory corrections")
        resps = pbm.normalize_phase(resps, min_phase=min_phase, max_phase=max_phase)

        tip_coords_resp = pbm.transform_xyz(tip_coords.copy(), resps, model_params)
        dist_coords_resp = pbm.transform_xyz(dist_coords.copy(), resps, model_params)
        prox_coords_resp = pbm.transform_xyz(prox_coords.copy(), resps, model_params)

        tip_coords_ukf_resp = pbm.transform_xyz(
            tip_coords_ukf.copy(), resps, model_params
        )
        dist_coords_ukf_resp = pbm.transform_xyz(
            dist_coords_ukf.copy(), resps, model_params
        )
        prox_coords_ukf_resp = pbm.transform_xyz(
            prox_coords_ukf.copy(), resps, model_params
        )

        inds = numpy.logical_and(trigs >= trig_window[0], trigs <= trig_window[1])
        trigs = trigs[inds]

        tip_coords = tip_coords[inds]
        dist_coords = dist_coords[inds]
        prox_coords = prox_coords[inds]

        tip_coords_resp = tip_coords_resp[inds]
        dist_coords_resp = dist_coords_resp[inds]
        prox_coords_resp = prox_coords_resp[inds]

        tip_coords_ukf = tip_coords_ukf[inds]
        dist_coords_ukf = dist_coords_ukf[inds]
        prox_coords_ukf = prox_coords_ukf[inds]

        tip_coords_ukf_resp = tip_coords_ukf_resp[inds]
        dist_coords_ukf_resp = dist_coords_ukf_resp[inds]
        prox_coords_ukf_resp = prox_coords_ukf_resp[inds]

        name = (output_dirname + "/{:04d}-sax-gated.mp4").format(frame_index)
        print('Writing SAx movie file: "{}"'.format(name))
        p = Plotter(
            tip_coords,
            dist_coords,
            prox_coords,
            tip_coords_ukf_resp,
            dist_coords_ukf_resp,
            prox_coords_ukf_resp,
            trigs,
            sax_cine_data,
            label="without correction",
            label_ukf="with correction",
        )
        p.render_to_file(name, fps=12)

        name = (output_dirname + "/{:04d}-lax-gated.mp4").format(frame_index)
        print('Writing LAx movie file: "{}"'.format(name))
        p = Plotter(
            tip_coords,
            dist_coords,
            prox_coords,
            tip_coords_ukf_resp,
            dist_coords_ukf_resp,
            prox_coords_ukf_resp,
            trigs,
            lax_cine_data,
            label="without correction",
            label_ukf="with correction",
        )
        p.render_to_file(name, fps=12)

        # name = (output_dirname + "/{:04d}-sax.mp4").format(frame_index)
        # print('Writing SAx movie file: "{}"'.format(name))
        # p = Plotter(
        #     tip_coords_resp,
        #     dist_coords_resp,
        #     prox_coords_resp,
        #     tip_coords_ukf_resp,
        #     dist_coords_ukf_resp,
        #     prox_coords_ukf_resp,
        #     trigs,
        #     sax_cine_data,
        # )
        # p.render_to_file(name, fps=12)

        # name = (output_dirname + "/{:04d}-lax.mp4").format(frame_index)
        # print('Writing LAx movie file: "{}"'.format(name))
        # p = Plotter(
        #     tip_coords_resp,
        #     dist_coords_resp,
        #     prox_coords_resp,
        #     tip_coords_ukf_resp,
        #     dist_coords_ukf_resp,
        #     prox_coords_ukf_resp,
        #     trigs,
        #     lax_cine_data,
        # )
        # p.render_to_file(name, fps=12)

        # name = (output_dirname + "/{:04d}-movie.mp4").format(frame_index)
        # print('Writing movie file: "{}"'.format(name))
        # p = Plotter4x4(
        #     (tip_coords, dist_coords, prox_coords),
        #     (tip_coords_ukf, dist_coords_ukf, prox_coords_ukf),
        #     (tip_coords_resp, dist_coords_resp, prox_coords_resp),
        #     (tip_coords_ukf_resp, dist_coords_ukf_resp, prox_coords_ukf_resp),
        #     trigs,
        #     sax_cine_data,
        #     lax_cine_data
        # )
        # p.render_to_file(name, fps=12)

        plt.close(fig="all")

    except KeyboardInterrupt:
        raise

    except:
        print("Something went wrong ... moving to the next frame")
        traceback.print_exc()
        continue
