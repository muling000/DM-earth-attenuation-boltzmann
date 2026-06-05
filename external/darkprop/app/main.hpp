#include <ios>
#include <map>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <csignal>
#include <mpi.h>
#include <spdlog/spdlog.h>
#include "darkprop/core.hpp"
#include "darkprop/Sun.hpp"
#include "darkprop/SolarDM.hpp"
#include "darkprop/SIHelmDM.hpp"
#include "darkprop/DMElectron.hpp"
#include "darkprop/DMHalo.hpp"
#include "darkprop/Injector.hpp"
#include "external/indicators.hpp"

using namespace darkprop;

enum class MediumModel { homoearth, homoelectronearth, prem, sun };
enum class DMModel {
    si, sihelm, solardm, dme,
#if __WITH_GENIE__
    geniebdm
#endif
};

inline std::string str_tolower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

inline DMModel parse_dm_model(const std::string& dm_model)
{
    std::string name = str_tolower(dm_model);
    if (name == "si" || name == "simple") { return DMModel::si; }
    if (name == "sihelm" || name == "elastic") { return DMModel::sihelm; }
    if (name == "solardm") { return DMModel::solardm; }
    if (name == "dmelectron") { return DMModel::dme; }
#if __WITH_GENIE__
    if (name == "geniebdm") { return DMModel::geniebdm; }
#endif
    throw std::invalid_argument(std::string("dm model ") + dm_model
                    + " unknown.\nsupported are \"si\" and \"sihelm\"");
}

inline MediumModel parse_medium_model(const std::string& medium_model)
{
    std::string name = str_tolower(medium_model);
    if (name == "homoearth" || name == "homo") { return MediumModel::homoearth; }
    if (name == "homoelectronearth") { return MediumModel::homoelectronearth; }
    if (name == "sun") { return MediumModel::sun; }
    if (name == "prem") { return MediumModel::prem; }
    throw std::invalid_argument(std::string("medium model ") + medium_model
    + " unknown.\nsupported are \"homoearth\", \"prem\", \"homoelectronearth\", and sun");
}


template<typename T>
std::vector<std::vector<std::vector<T>>> make_triple_vector(
        std::size_t dim1, std::size_t dim2, std::size_t dim3=0)
{
    std::vector<std::vector<std::vector<T>>> v3d(dim1, std::vector<std::vector<T>>(dim2));
    if (dim3 > 0) {
        for (auto& v2d : v3d) {
            std::fill(v2d.begin(), v2d.end(), std::vector<T>(dim3));
        }
    }
    return v3d;
}

template<typename Vector3, typename Value, typename TStore>
std::size_t ring_classify(
        std::vector<Event<Vector3, Value>>& events,
        const std::vector<Value>& ring_bins,
        const Vector3& axis,
        std::vector<std::vector<Event<Vector3, Value>>>& ring_events,
        std::vector<std::vector<TStore>>& ring_weights,
        const std::vector<std::size_t>& sample_nums,
        std::size_t sample_sz)
{
    std::size_t sample_num = 0;
    for (auto& event : events) {
        Value theta = std::acos(event.r.dot(axis) / event.r.norm());
        if (theta <= ring_bins.front() || theta > ring_bins.back()) { continue; }
        std::size_t ring = std::distance(
            ring_bins.begin() + 1,
            std::lower_bound(ring_bins.begin() + 1, ring_bins.end(), theta));
        if (sample_nums[ring] + ring_events[ring].size() < sample_sz) {
            ring_events[ring].push_back(event);
            ring_weights[ring].push_back(static_cast<TStore>(event.weight));
            ++sample_num;
        }
    }
    return sample_num;
}

inline
bool continue_run(const std::vector<std::vector<std::size_t>>& sample_nums,
                  std::size_t sample_sz) {
    for (const auto& d1 : sample_nums) {
        for (const auto& size : d1) {
            if (size < sample_sz) return true;
        }
    }
    return false;
}

template<typename Value>
std::string make_filename(const std::string& dir, const std::string& id,
                          Value mchi, Value sigma,
                          const std::string& what="event",
                          std::size_t track=0)
{
    std::stringstream ss;
    ss.setf(std::ios_base::scientific);
    ss.precision(3);
    if (what == "event") {
        ss << dir << "/" << id << "_mchi" << mchi
            << "GeV_sigma" << sigma / units::cm / units::cm << "cm2.hdf5";
    } else if (what == "track") {
        ss << dir << "/" << id << "_mchi" << mchi
            << "GeV_sigma" << sigma / units::cm / units::cm
            << "cm2_fulltrack_" << track << ".csv";
    }
    return ss.str();
}

template<> inline
std::string make_filename<std::size_t>(const std::string& dir,
                                       const std::string& id,
                                       std::size_t mi, std::size_t si,
                                       const std::string& what,
                                       std::size_t track)
{
    std::stringstream ss;
    if (what == "event") {
        ss << dir << "/" << id << "_m" << mi << "_s" << si << ".hdf5";
    } else if (what == "track") {
        ss << dir << "/" << id << "_m" << mi << "_s" << si
            << "_fulltrack_" << track << ".csv";
    }
    return ss.str();
}


inline std::size_t my_nums(int my_rank, int comm_sz, std::size_t total)
{
    std::size_t base_num {total / static_cast<std::size_t>(comm_sz)};
    std::size_t remaind {total % static_cast<std::size_t>(comm_sz)};
    return base_num + (static_cast<std::size_t>(my_rank) < remaind ? 1 : 0);
}

inline std::size_t my_offset(int my_rank, int comm_sz, std::size_t total)
{
    std::size_t base_num {total / static_cast<std::size_t>(comm_sz)};
    std::size_t remaind {total % static_cast<std::size_t>(comm_sz)};
    return static_cast<std::size_t>(my_rank) < remaind ?
        my_rank * (base_num + 1) : remaind + my_rank * base_num;
}

inline std::size_t get_width(std::size_t size)
{
    std::size_t w {2};
    while (static_cast<std::size_t>(std::pow(10, w - 1)) / (size + 1) == 0) { ++w; }
    return w;
}

template<typename TStore=float, typename Vector3, typename Value>
std::vector<std::vector<std::vector<std::size_t>>> run(
    std::map<std::string, hid_t> dsets,
    std::shared_ptr<Particle<Vector3, Value>> dmp,
    std::shared_ptr<MediumBall<Vector3, Value>> mp,
    std::shared_ptr<Injector<Vector3, Value>> injp,
    std::vector<Value> depth,
    std::vector<std::string> depth_groups,
    std::vector<std::string> ring_groups,
    Value Tcut,
    int t_rank,
    int t_comm_sz,
    std::size_t total_samples,
    std::size_t maximal_simulation,
    std::size_t full_tracks,
    const std::string& track_dir,
    const std::string& id,
    Vector3 axis,
    const std::vector<Value>& ring_bins,
    bool store_time,
    int seed,
    DMModel dmmodel,
    std::size_t chunk_size,
    std::size_t loglag,
    unsigned long long_track)
{
    std::stringstream srank;
    srank << std::setw(get_width(t_comm_sz)) << t_rank;
    indicators::ProgressBar bar{
        indicators::option::BarWidth{40},
        indicators::option::Start{" ["},
        indicators::option::Fill{"█"},
        indicators::option::Lead{"█"},
        indicators::option::Remainder{"-"},
        indicators::option::End{"]"},
        indicators::option::ShowPercentage(true),
        indicators::option::PrefixText{"Process" + srank.str()},
        indicators::option::ForegroundColor{indicators::Color::white},
        indicators::option::ShowElapsedTime{true},
        indicators::option::FontStyles{
            std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
    };
    RandomNumber<Value> rn(seed);

    const std::size_t tracks = my_nums(t_rank, t_comm_sz, full_tracks);
    const std::size_t track_offset = my_offset(t_rank, t_comm_sz, full_tracks);
    std::size_t track = 0;
    while (track < tracks) {
        std::string track_f = make_filename(
            track_dir, id, dmp->m, dmp->sigma0, "track", track + track_offset);
        injp->next(*dmp);
        dump_track_csv(simulate_track(*dmp, *mp, Tcut, rn), track_f);
        ++track;
    }

    std::size_t rings = ring_bins.size() - 1;
    auto cross_events = make_triple_vector<Event<Vector3, Value>>(depth.size(), rings);
    auto weights = make_triple_vector<TStore>(depth.size(), rings);

    std::vector<std::vector<std::size_t>> sample_nums(depth.size());
    std::fill(sample_nums.begin(), sample_nums.end(), std::vector<std::size_t>(rings, 0));
    // {Nsample, Nsim}
    auto counts = make_triple_vector<std::size_t>(depth.size(), rings, 2);

    const std::size_t sample_sz = my_nums(t_rank, t_comm_sz, total_samples);
    const std::size_t baseoffset = my_offset(t_rank, t_comm_sz, total_samples);
    const std::size_t max_runnum = my_nums(t_rank, t_comm_sz, maximal_simulation);

    std::size_t runnum = 0;
    std::size_t total_sample_num = 0;
    bool keep_going = true;

    signal(SIGINT, [](int signum) {
        indicators::show_console_cursor(true);
        std::cout << termcolor::reset << std::endl;
        exit(signum);
    });
    indicators::show_console_cursor(false);

    do {
        auto weight = injp->next(*dmp);

        if (dmmodel == DMModel::si) { mp->setCache(*dmp); }

        auto events = simulate_cross_depth(
            *dmp, *mp, Tcut, depth, rn, weight, long_track);
        ++runnum;

        for (std::size_t di = 0; di < depth.size(); ++di) {
            total_sample_num += ring_classify(events[di],
                                              ring_bins,
                                              axis.normalized(),
                                              cross_events[di],
                                              weights[di],
                                              sample_nums[di],
                                              sample_sz);
            for (std::size_t ri = 0; ri != rings; ++ri) {
                // set chunk size for balance between memory usage and disk IO
                std::size_t events_on_ring = cross_events[di][ri].size();
                if (sample_nums[di][ri] == sample_sz) { continue; }
                if (sample_nums[di][ri] > sample_sz) {
                    throw std::runtime_error("sample_nums > sample_sz");
                }
                if (sample_nums[di][ri] + events_on_ring == sample_sz
                    || events_on_ring > chunk_size || runnum == max_runnum) {
                    auto groupname = rings > 1 ?
                                     depth_groups[di] + "/" + ring_groups[ri]
                                     : depth_groups[di];
                    std::size_t offset = baseoffset + sample_nums[di][ri];
                    h5dump_events_galactic<TStore>(dsets, groupname,
                        cross_events[di][ri], offset, store_time);
                    write_dset(dsets[groupname + "/weight"], weights[di][ri],
                               {offset}, {weights[di][ri].size()});

                    sample_nums[di][ri] += events_on_ring;

                    if (sample_nums[di][ri] == sample_sz
                        || runnum == max_runnum) {
                        counts[di][ri][0] = sample_nums[di][ri];
                        counts[di][ri][1] = runnum;
                    }

                    cross_events[di][ri].clear();
                    weights[di][ri].clear();
                }
            }
        }

        keep_going = continue_run(sample_nums, sample_sz) && runnum < max_runnum;
        if (runnum % loglag == 0 || !keep_going) {
            bar.set_option(indicators::option::PostfixText{"run # " + std::to_string(runnum)});
            auto progress = total_sample_num * 100.0 / (sample_sz * rings * depth.size());
            bar.set_progress(progress);
        }

    } while (keep_going);
    indicators::show_console_cursor(true);
    return counts;
}

template<typename Value>
void load_intensity(const std::string& filename,
                    const std::string& mgroup,
                    std::vector<Value>& t_b,
                    std::vector<Value>& t_l,
                    std::vector<Value>& t_T,
                    std::vector<Value>& t_I)
{
    hid_t fd = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    t_b = read_dset<Value>(fd, mgroup + "/b");
    t_l = read_dset<Value>(fd, mgroup + "/l");
    t_T = read_dset<Value>(fd, mgroup + "/T");
    t_I = read_dset<Value>(fd, mgroup + "/intensity");
    H5Fclose(fd);
}
