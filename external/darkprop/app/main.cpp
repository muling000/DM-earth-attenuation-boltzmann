#include <chrono>
#include <memory>
#include <sstream>
#include <iostream>
#include <filesystem>

#include "main.hpp"
#include "external/yuc_config.h"

#if __WITH_GENIE__
#include "darkprop/GENIEBDM.hpp"
#endif

using Value = double;
using TStore = float;
using Vector3 = Vector3d<Value>;


int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cout << "Usage: mpiexec -n N " << argv[0] << " <config.toml>" << std::endl;
        exit(1);
    }

    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    {

    yuc::config config;
    config.parse_toml(argv[1]);

    const std::string id = config["id"];

    Value Tmin {0}, Tmax {0};
    // [input]
    const auto input = config["input"];
    const std::string input_type = input["type"];
    std::vector<Value> mchis = input["mchis"];
    std::vector<Value> sigmas = input["sigmas"];
    for (auto& s : sigmas) { s *= units::cm*units::cm; }
    const std::string mchi_file = input["mchi_file"].or_get(std::string(""));
    const std::string sigma_file = input["sigma_file"].or_get(std::string(""));
    const Value mmed = input["mmed"].or_get(0.0);
    const std::string mediator_limit = input["mediator_limit"].or_get(std::string(""));
    std::ifstream ifs;
    if (std::filesystem::exists(mchi_file)) {
        ifs.open(mchi_file);
        // overwrite mchis
        Value tmp;
        mchis.clear();
        while (ifs >> tmp) { mchis.push_back(tmp); }
        ifs.close();
    }
    if (std::filesystem::exists(sigma_file)) {
        ifs.open(sigma_file);
        // overwrite sigmas
        Value tmp;
        sigmas.clear();
        while (ifs >> tmp) { sigmas.push_back(tmp * units::cm * units::cm); }
        ifs.close();
    }
    Value Tcut = input["Tcut"].or_get(0.0);

    // [input.halo]
    const auto input_halo = input["halo"];
    std::vector<Value> vec_earth {0, 0, 0};
    // for g++ < 10
    if (input_halo["v_earth"].is_set()) { vec_earth = input_halo["v_earth"]; }

    const Value kmsec = units::km / units::sec;
    const Vector3 v_earth(vec_earth[0]*kmsec, vec_earth[1]*kmsec, vec_earth[2]*kmsec);
    const Value halo_v_esc = input_halo["v_esc"].num(544) * kmsec;
    const Value halo_v_0 = input_halo["v_0"].num(220) * kmsec;
    const Value halo_v_cut = input_halo["v_cut"].num(1e-5) * kmsec;
    const Value halo_v_min = input_halo["v_min"].num(1e-5) * kmsec;
    Value halo_v_max = input_halo["v_max"].num(-1) * kmsec;
    halo_v_max = halo_v_max > 0 ? halo_v_max : halo_v_esc + v_earth.norm();

    // [input.flux]
    auto input_flux = input["flux"];
    bool point_source = input_flux["point_source"].or_get(false);
    const std::string flux_file = input_flux["flux_file"].or_get(std::string(""));
    std::vector<Value> energy;
    std::vector<std::vector<Value>> fluxes(mchis.size()*sigmas.size());
    if (input_type == "flux") {
        Value e_tmp;
        Value flux_tmp;
        // TODO check the file
        ifs.open(flux_file);
        while (ifs >> e_tmp) {
            energy.push_back(e_tmp);
            for (std::size_t i = 0; i != fluxes.size(); ++i) {
                ifs >> flux_tmp;
                fluxes[i].push_back(flux_tmp);
            }
        }
        ifs.close();

        Tmin = input_flux["Tmin"].or_get(0.0);
        Tmax = input_flux["Tmax"].or_get(0.0);
        if (Tmin <= 0.0) {
            Tmin = *std::min_element(energy.cbegin(), energy.cend());
        }
        if (Tmax <= 0.0) {
            Tmax = *std::max_element(energy.cbegin(), energy.cend());
        }
        if (Tcut <= 0.0) {
            Tcut = Tmin;
        }
    }

    // [input.intensity]
    auto input_intensity = input["intensity"];
    const std::string intensity_file
        = input_intensity["intensity_file"].or_get(std::string(""));
    std::vector<std::size_t> select_mass = input_intensity["select_mass"].or_get({});
    Value base_sigma {0};
    if (input_type == "intensity") {
        mchis = read_dset<Value>(intensity_file, "mchi");
        if (select_mass.size() != 0) {
            auto m_tmp = mchis;
            mchis.resize(select_mass.size());
            for (std::size_t i = 0; i != select_mass.size(); ++i) {
                if (select_mass[i] >= m_tmp.size()) {
                    throw std::out_of_range("select_mass index out of range");
                }
                mchis[i] = m_tmp[select_mass[i]];
            }
        } else {
            for (std::size_t i = 0; i != mchis.size(); ++i) {
                select_mass.push_back(i);
            }
        }
        base_sigma = read_dset<Value>(intensity_file, "sigma")[0];
    }

    // [input.sample]
    auto input_sample = input["sample"];
    const std::string sample_file = input_sample["sample_file"].or_get(std::string(""));
    std::size_t my_inj_offset {0};
    std::size_t my_inj_num {0};
    if (input_type == "sample") {
        hid_t fd = H5Fopen(sample_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t ds = H5Dopen(fd, "ran_E", H5P_DEFAULT);
        hsize_t dims[2];
        hid_t dspace = H5Dget_space(ds);
        H5Sget_simple_extent_dims(dspace, dims, NULL);
        std::size_t input_sample_size = static_cast<std::size_t>(dims[0]);

        my_inj_num = my_nums(my_rank, comm_sz, input_sample_size);
        my_inj_offset = my_offset(my_rank, comm_sz, input_sample_size);

        mchis = read_dset<Value>(fd, "mchi");
        H5Sclose(dspace);
        H5Dclose(ds);
        H5Fclose(fd);
    }

    // [input.source]
    auto input_source = input["source"];
    const std::string source_inv_cdf_r_file
        = input_source["inv_cdf_r_file"].or_get(std::string(""));
    const std::string source_inv_cdf_T_file
        = input_source["inv_cdf_T_file"].or_get(std::string(""));

    // [input.solardm]
    auto input_solardm = input["solardm"];
    const std::string solar_density_integral_file
        = input_solardm["solar_density_integral_file"].or_get(std::string(""));
    Value temperature_scale = input_solardm["temperature_scale"].or_get(1.0);

#if __WITH_GENIE__
    // [input.geniebdm]
    auto input_geniebdm = input["geniebdm"];
    const std::string genie_tune = input_geniebdm["tune"].or_get(std::string(""));
    const std::string genie_event_gen_list
        = input_geniebdm["event_generator_list"].or_get(std::string(""));
    const std::vector<std::string> genie_xsec_files
        = input_geniebdm["xsec_files"].or_get({std::string("")});
    if (genie_xsec_files.front() != "" && mchis.size() != genie_xsec_files.size()) {
        throw std::invalid_argument("number of xsec files != number of DM masses");
    }
    const std::string genie_cache_file
        = input_geniebdm["cache_file"].or_get(std::string(""));
    const std::string genie_mesg_thrs
        = input_geniebdm["messenger_thresholds"].or_get(std::string(""));
    const Value origin_gz = input_geniebdm["origin_gz"].or_get(1.0);
    const std::size_t xsec_spline_knots
        = input_geniebdm["xsec_spline_knots"].or_get(1000L);
#endif

    // [output]
    const auto output = config["output"];
    const std::size_t sample_size = output["sample_size"];
    const std::string outdir = output["outdir"];
    const bool overwrite_output = output["overwrite"].or_get(false);
    const std::string file_name_style = output["file_name_style"].or_get(std::string("long"));
    std::size_t rings = output["rings"].or_get(1L);
    rings = rings < 1 ? 1 : rings;
    const Value cos_theta_lower = output["cos_theta_lower"].or_get(-1.0);
    const Value cos_theta_upper = output["cos_theta_upper"].or_get(1.0);
    if (cos_theta_lower >= cos_theta_upper) {
        throw std::invalid_argument("cos_theta_lower >= cos_theta_upper");
    }
    const Value theta_lower = std::acos(cos_theta_upper);
    const Value theta_upper = std::acos(cos_theta_lower);
    const std::string ring_space = output["ring_space"].or_get(std::string("linear"));
    const std::size_t chunk_size = output["chunk_size"].or_get(10000L);
    const std::vector<Value> r_axis = output["axis"].or_get({1.0, 0.0, 0.0});
    const bool store_time = output["store_time"].or_get(false);

    Vector3 ring_axis;
    if (input_type == "halo") {
        ring_axis = {v_earth[0], v_earth[1], v_earth[2]};
    } else {
        ring_axis = {r_axis[0], r_axis[1], r_axis[2]};
    }
    const std::size_t maximal_simulation = output["maximal_simulation"].or_get(1000000000L);

    const std::size_t full_tracks = output["full_tracks"].or_get(0L);

    // [output.log]
    auto output_log = output["log"];
    std::size_t screen_log_lag = output_log["screen_log_lag"].or_get(10000L);
    unsigned long warning_long_track = output_log["warning_long_track"].or_get(1000000L);

    // [numerical]
    auto numerical = config["numerical"];
    int rn_seed = numerical["random_seed"].or_get(-1L);
    rn_seed = rn_seed >= 0 ? rn_seed + my_rank : rn_seed;
    const Value interp_Tmin = numerical["interp_Tmin"].or_get(1e-100);
    const std::size_t interp_num = numerical["interp_num"].or_get(10000L);
    const std::size_t sample_buffer_size
        = numerical["sample_buffer_size"].or_get(10000L);

    // initialization
    // depth
    std::vector<Value> depth = output["depth"];

    // depth group name
    std::stringstream ss;
    std::vector<std::string> depth_groups;
    for (std::size_t di = 0; di != depth.size(); ++di) {
        ss << "depth_" << di;
        depth_groups.push_back(ss.str());
        depth[di] *= units::km;
        ss.str("");
    }

    // ring group name
    std::vector<std::string> ring_groups;
    for (std::size_t r = 0; r != rings; ++r) {
        ss << "ring_" << r;
        ring_groups.push_back(ss.str());
        ss.str("");
    }

    std::vector<Value> ring_bins;
    if (ring_space == "linear") {
        ring_bins = linspace<Value>(theta_lower, theta_upper, rings + 1, true);
    } else if (ring_space == "cos_linear") {
        ring_bins = linspace<Value>(cos_theta_upper, cos_theta_lower,
                                    rings + 1, true);
        for (auto& r : ring_bins) { r = std::acos(r); }
    } else {
        throw std::invalid_argument(
            std::string("ring_space value ")
            + ring_space
            + " unknown.\nsupported are \"linear\" and \"cos_linear\"");
    }

    // DM
    std::shared_ptr<Particle<Vector3, Value>> dmp;
    const std::string dm_model = input["dm_model"];
    DMModel dmmodel = parse_dm_model(dm_model);

    // Medium
    std::shared_ptr<MediumBall<Vector3, Value>> mp;
    auto premp = std::make_shared<PREMEarth<Vector3, Value>>();
    const std::string medium_model = input["medium_model"];
    MediumModel mediummodel = parse_medium_model(medium_model);
    switch (mediummodel) {
        case MediumModel::homoearth:
            mp = std::make_shared<HomoEarth<Vector3, Value>>();
            break;
        case MediumModel::prem:
            mp = premp;
            break;
        case MediumModel::homoelectronearth:
            mp = std::make_shared<HomoElectronEarth<Vector3, Value>>();
            break;
        case MediumModel::sun:
            auto sunp = std::make_shared<Sun<Vector3, Value>>();
            sunp->init(solar_density_integral_file);
            mp = sunp;
            break;
    }

    // injector
    std::shared_ptr<Injector<Vector3, Value>> injp;
    auto flux_injp = std::make_shared<FluxInjector<Vector3, Value>>(rn_seed);
    auto sample_injp
        = std::make_shared<SampleInjector<Vector3, Value>>(rn_seed);
    auto intensity_injp
        = std::make_shared<IntensityInjector<Vector3, Value>>(rn_seed);
    if (input_type == "halo") {
        auto halo_injp = std::make_shared<HaloInjector<Vector3, Value>>(
            halo_v_min, halo_v_max, v_earth, halo_v_esc, halo_v_0,
            mp->getRadius(), 0.0, rn_seed);
        injp = halo_injp;
    } else if (input_type == "flux") {
        injp = flux_injp;
    } else if (input_type == "intensity") {
        intensity_injp->setBaseSigma(base_sigma);
        injp = intensity_injp;
    } else if (input_type == "sample") {
        injp = sample_injp;
    } else if (input_type == "source") {
        injp = std::make_shared<SourceInjector<Vector3, Value>>(
            source_inv_cdf_r_file,
            source_inv_cdf_T_file,
            "linear",
            rn_seed);
    } else {
        throw std::invalid_argument(
            std::string("input.type value ")
            + input_type
            + " unknown.\nsupported are \"halo\","
            + "\"flux\", \"intensity\", \"sample\", and \"source\"."
        );
    }

    auto counts = make_triple_vector<std::size_t>(depth.size(), rings, 2);

    std::string dir = outdir + "/" + id;
    std::filesystem::create_directories(dir);
    std::string track_dir = dir + "/fulltracks";
    if (full_tracks > 0) { std::filesystem::create_directories(track_dir); }

    auto start_time = std::chrono::high_resolution_clock::now();
    for (std::size_t mi = 0; mi != mchis.size(); mi++) {
        Value mchi = mchis[mi];
        if (input_type == "intensity") {
            std::vector<Value> b, l, T, I;
            load_intensity(intensity_file,
                           std::string("m_")+std::to_string(select_mass[mi]),
                           b, l, T, I);

            intensity_injp->init(b, l, T, I);
            Tmin = *std::min_element(T.cbegin(), T.cend());
            Tmax = *std::max_element(T.cbegin(), T.cend());
            // TODO user set Tcut
            Tcut = Tmin;
        }
        for (std::size_t si = 0; si != sigmas.size(); ++si) {
            Value sigma = sigmas[si];
            if (input_type == "halo") {
                Tmax = 0.5 * mchi * halo_v_max * halo_v_max;
                Tmin = 0.5 * mchi * halo_v_min * halo_v_min;
                Tcut = 0.5 * mchi * halo_v_cut * halo_v_cut;
            } else if (input_type == "flux") {
                flux_injp->init(energy, fluxes[si + mi * sigmas.size()]);
                flux_injp->m_Tmin = Tmin;
                flux_injp->m_Tmax = Tmax;
                flux_injp->m_point_source = point_source;
            } else if (input_type == "intensity") {
                intensity_injp->setSigma(sigma/units::cm/units::cm);
            } else if (input_type == "sample") {
                sample_injp->open(sample_file, H5F_ACC_RDONLY);
                sample_injp->setBufferSize(sample_buffer_size);
                sample_injp->setWeightMassIdx(mi);
                sample_injp->setDataRange(my_inj_offset, my_inj_offset + my_inj_num);
            }

            switch (dmmodel) {
                case DMModel::si:
                    dmp = std::make_shared<SIDM<Vector3, Value>>(sigma, mchi);
                    break;
                case DMModel::sihelm:
                {
                    auto sihelmdmp =
                        std::make_shared<SIHelmDM<Vector3, Value>>(sigma, mchi);
                    for (auto& t : mp->getTargets()) {
                        sihelmdmp->initCrossSectionCDF(interp_Tmin, Tmax,
                                                        t, interp_num);
                    }
                    premp->cross_section_energy_dependent = true;
                    dmp = sihelmdmp;
                    break;
                }
                case DMModel::solardm:
                {
                    auto solardmp = std::make_shared<SolarDM<Vector3, Value>>(sigma, mchi);
                    solardmp->setTempScale(temperature_scale);
                    dmp = solardmp;
                    break;
                }
                case DMModel::dme:
                {
                    auto dmep = std::make_shared<DMElectron<Vector3, Value>>(
                                    sigma, mchi, mmed);
                    dmep->mediator_limit = mediator_limit;
                    dmp = dmep;
                    break;
                }
#if __WITH_GENIE__
                case DMModel::geniebdm:
                {
                    genie_init(genie_tune, genie_mesg_thrs, genie_xsec_files[mi],
                               genie_event_gen_list, rn_seed, genie_cache_file);
                    genie_add_dark_matter(mchi, mmed);
                    genie_set_coupling(origin_gz);
                    auto geniebdmp = std::make_shared<GENIEBDM<Vector3, Value>>(
                        sigma, 0.0, Tmin, Tmax, origin_gz);
                    geniebdmp->initGEVGDriver(*mp, xsec_spline_knots);
                    dmp = geniebdmp;
                    break;
                }
#endif
            }

            MPI_Barrier(MPI_COMM_WORLD);
            auto filename = make_filename(dir, id, mchi, sigma);
            if (file_name_style == "short") {
                if (input_type == "intensity") {
                    filename = make_filename(dir, id, select_mass[mi], si);
                } else {
                    filename = make_filename(dir, id, mi, si);
                }
            }
            if (std::filesystem::exists(filename)) {
                if (overwrite_output) {
                    if (my_rank == 0) { spdlog::warn("overwriting {}", filename); }
                } else {
                    if (my_rank == 0) { spdlog::info("{} exists, skipped", filename); }
                    continue;
                }
            }
            if (my_rank == 0) {
                spdlog::info("Collectively generating file: {} ....", filename);
            }
            hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
            auto fd = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
            H5Pclose(fapl);

            hid_t ds = make_dset<TStore>(fd, "sigma");
            write_dset(ds, sigma / units::cm / units::cm);
            write_attr(ds, "unit", "cm^2");
            H5Dclose(ds);

            ds = make_dset<TStore>(fd, "mchi");
            write_dset(ds, mchi);
            write_attr(ds, "unit", "GeV");
            H5Dclose(ds);

            ds = make_dset<TStore>(fd, "Tcut");
            write_dset(ds, Tcut);
            write_attr(ds, "unit", "GeV");
            H5Dclose(ds);
            if (input_type != "sample") {
                ds = make_dset<TStore>(fd, "Tmin");
                write_dset(ds, Tmin);
                write_attr(ds, "unit", "GeV");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "Tmax");
                write_dset(ds, Tmax);
                write_attr(ds, "unit", "GeV");
                H5Dclose(ds);
            }
            ds = make_dset<TStore>(fd, "R0");
            write_dset(ds, mp->getRadius() / units::km);
            write_attr(ds, "unit", "km");
            H5Dclose(ds);
            if (input_type == "halo") {
                ds = make_dset<TStore>(fd, "Vesc");
                write_dset(ds, halo_v_esc / kmsec);
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "V0");
                write_dset(ds, halo_v_0 / kmsec);
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "Vearth", {vec_earth.size()});
                write_dset(ds, vec_earth, {0}, {vec_earth.size()});
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "Vmin");
                write_dset(ds, halo_v_min / kmsec);
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "Vmax");
                write_dset(ds, halo_v_max / kmsec);
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
                ds = make_dset<TStore>(fd, "Vcut");
                write_dset(ds, halo_v_cut / kmsec);
                write_attr(ds, "unit", "km/s");
                H5Dclose(ds);
            } else if (input_type == "flux") {
                ds = make_dset<TStore>(fd, "norm");
                write_dset(ds, flux_injp->m_norm);
                write_attr(ds, "unit", "cm^-2 sr^-1");
                H5Dclose(ds);
            } else if (input_type == "intensity") {
                ds = make_dset<TStore>(fd, "norm");
                write_dset(ds, intensity_injp->getNorm());
                write_attr(ds, "unit", "cm^-2 sr^-1");
                H5Dclose(ds);
            }
            auto dsets = make_dsets(
                fd, depth_groups, ring_groups, depth, ring_bins,
                static_cast<Value>(mp->getRadius()), sample_size, chunk_size, store_time);
            if (my_rank == 0) {
                spdlog::info("Collectively generating file: {} done\nStart simulation",
                             filename);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            auto my_counts = run(dsets, dmp, mp, injp, depth, depth_groups,
                                 ring_groups, Tcut, my_rank, comm_sz,
                                 sample_size, maximal_simulation, full_tracks,
                                 track_dir, id, ring_axis, ring_bins,
                                 store_time, rn_seed, dmmodel, chunk_size,
                                 screen_log_lag, warning_long_track);
            for (std::size_t di = 0; di != depth.size(); ++di) {
                for (std::size_t ri = 0; ri != rings; ++ri) {
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Reduce(my_counts[di][ri].data(), counts[di][ri].data(),
                               2, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                               0, MPI_COMM_WORLD);
                    auto groupname = rings > 1 ?
                                     depth_groups[di] + "/" + ring_groups[ri]
                                     : depth_groups[di];
                    if (my_rank == 0) {
                        write_dset(dsets[groupname + "/Nsample"], counts[di][ri][0]);
                        write_dset(dsets[groupname + "/Nsim"], counts[di][ri][1]);
                    }
                }
            }
            close_dsets(dsets);
            H5Fclose(fd);
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();

    if (my_rank == 0) {
        spdlog::info("costs {} seconds", std::chrono::duration_cast<std::chrono::seconds>(
                end_time - start_time).count());
    }

    }  // release resources before MPI_Finalize();
    MPI_Finalize();
    return 0;
}
