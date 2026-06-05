/**
 * @file
 * Injectors for injecting dark matter particle.
 */
#pragma once
#include <functional>
#include "Base.hpp"
#include "IO.hpp"
#include "Vector3d.hpp"
#include "numerical/GSLInterpolator.hpp"
#include "numerical/GridLinearInterpolator.hpp"
#include "numerical/Integrator.hpp"

namespace darkprop {

using namespace numerical;

/**
 * The abstract base class for injecting DM.
 */
template<typename Vector3, typename Value>
class Injector {
public:
    Value m_radius;
    Value m_t0;
    RandomNumber<Value> m_rn;
    virtual ~Injector() {}
    virtual Value next(Particle<Vector3, Value>&) = 0;
    virtual Value operator()(Particle<Vector3, Value>& p) {
        return next(p);
    }
    void setSeed(int t_seed) { m_rn = RandomNumber<Value>(t_seed); }
    void sett0(Value t_0) { m_t0 = t_0; }
};

template<typename Vector3, typename Value>
class HaloInjector : public Injector<Vector3, Value>
{
public:
    Value m_vmin;
    Value m_vmax;
    Vector3 m_vearth;
    Value m_vesc;
    Value m_v0;
    Value m_norm;
    explicit HaloInjector(Value t_vmin=0,
                          Value t_vmax=-1,
                          Vector3 t_vearth={0, 0, 240*units::km/units::sec},
                          Value t_vesc=constants::vesc,
                          Value t_v0=constants::v0,
                          Value t_radius=constants::rEarth,
                          Value t_t0=0,
                          int t_seed=-1)
        : m_vmin {t_vmin}, m_vmax {t_vmax}, m_vearth {t_vearth},
          m_vesc {t_vesc}, m_v0 {t_v0}
        {
            this->m_radius = t_radius;
            this->m_t0 = t_t0;
            this->m_rn = RandomNumber<Value>(t_seed);
            set_norm();
        }
    void set_norm() {
        Value vmax = m_vmax > 0 ? m_vmax : m_vearth.norm() + m_vesc;
        Integrator<std::function<Value(Value)>> integrator(true);
        m_norm = integrator([this](Value v) -> Value {
            return fv_halo_1d_earth(v, this->m_vearth.norm(), this->m_vesc,
                                    this->m_v0); }, m_vmin, vmax);
    }
    Value next(Particle<Vector3, Value>& p) override {
        init_halo(p, this->m_t0, m_vmin, m_vearth, this->m_rn,
                  1.1*this->m_radius, m_vmax, m_vesc, m_v0, this->m_radius);
        inject(p, this->m_radius);
        return p.v.norm() * m_norm;
    }
};

template<typename Vector3, typename Value>
class FluxInjector : public Injector<Vector3, Value>
{
public:
    Value m_point_source;
    Value m_Tmin;
    Value m_Tmax;
    Value m_norm;  ///< [cm^-2 sr^-1]
    GSLInterpolator<Value> spectrum;
    explicit FluxInjector(int seed=-1, Value t_radius=constants::rEarth) {
        this->setSeed(seed);
        this->m_radius = t_radius;
    }
    FluxInjector(const std::vector<Value>& t_T,
                 const std::vector<Value>& t_flux,
                 Value t_radius=constants::rEarth,
                 int seed=-1) {
        this->m_radius = t_radius;
        this->setSeed(seed);
        init(t_T, t_flux);
    }
    void init(const std::vector<Value>& t_T,
              const std::vector<Value>& t_flux) {
        m_Tmin = *std::min_element(t_T.begin(), t_T.end());
        m_Tmax = *std::max_element(t_T.begin(), t_T.end());
        spectrum.interpolate(t_T, t_flux);
        Integrator<GSLInterpolator<Value>> integrator(true);

        m_norm = integrator(spectrum, m_Tmin, m_Tmax);
    }
    ~FluxInjector() {}
    Value next(Particle<Vector3, Value>& p) override {
        Value T = std::exp(this->m_rn.uniform_ab(std::log(m_Tmin),
                                                 std::log(m_Tmax)));
        if (this->m_point_source) {
            Vector3 ep(0.0, 0.0, -1.0);
            init_T_point_source(
                p, T, ep, this->m_rn, this->m_t0, 1.1*this->m_radius, this->m_radius);
        } else {
            init_T_isotropic(
                p, T, this->m_rn, this->m_t0, 1.1*this->m_radius, this->m_radius);
        }
        inject(p, this->m_radius);
        return spectrum(T) * T * std::log(m_Tmax / m_Tmin) / m_norm;
    }
};

template<typename Vector3, typename Value>
class IntensityInjector : public Injector<Vector3, Value>
{
private:
    GridLinearInterpolator<3, std::vector<Value>> intensity_core;
    Value m_Tmin, m_Tmax, m_bmin, m_bmax, m_lmin, m_lmax;
    Value m_norm;   ///< [cm^-2 sr^-1]
    Value m_sigma;  ///< [cm^2]
    Value m_base_sigma;
public:
    explicit IntensityInjector(int seed=-1, Value t_radius=constants::rEarth)
        : m_Tmin {0}, m_Tmax {0}, m_bmin {0}, m_bmax{0},
          m_lmin {0}, m_lmax {0}, m_norm {0}, m_sigma {0}, m_base_sigma {0} {
        this->setSeed(seed);
        this->m_radius = t_radius;
    }
    IntensityInjector(const std::vector<Value>& t_b,
                      const std::vector<Value>& t_l,
                      const std::vector<Value>& t_T,
                      const std::vector<Value>& t_I,
                      int seed=-1, Value t_radius=constants::rEarth)
        : IntensityInjector(seed, t_radius) {
        init(t_b, t_l, t_T, t_I);
    }
    void init(const std::vector<Value>& t_b,
              const std::vector<Value>& t_l,
              const std::vector<Value>& t_T,
              const std::vector<Value>& t_I) {
        std::vector<Value> logT(t_T.size());
        std::vector<Value> logI(t_I.size());
        std::transform(t_T.cbegin(), t_T.cend(), logT.begin(),
            [](const auto& t) {return std::log(t);});
        std::transform(t_I.cbegin(), t_I.cend(), logI.begin(),
            [](const auto& i) {return i > 0 ? std::log(i) : std::log(1e-300);});
        intensity_core.interpolate({t_b, t_l, logT}, logI);
        m_bmin = intensity_core.getXmin(0);
        m_bmax = intensity_core.getXmax(0);
        m_lmin = intensity_core.getXmin(1);
        m_lmax = intensity_core.getXmax(1);
        m_Tmin = std::exp(intensity_core.getXmin(2));
        m_Tmax = std::exp(intensity_core.getXmax(2));

        MCIntegrator<std::function<Value(double *, std::size_t)>> integrator;

        m_norm = integrator([this](double *x, std::size_t dim) -> Value {
            (void) dim;
            return std::exp(this->intensity_core(x, dim) + x[2]) * std::cos(x[0]);},
            {m_bmin, m_lmin, std::log(m_Tmin)},
            {m_bmax, m_lmax, std::log(m_Tmax)});
    }
    void setSigma(Value t_sigma) { m_sigma = t_sigma; };
    void setBaseSigma(Value t_base_sigma) { m_base_sigma = t_base_sigma; };
    Value getNorm() { return m_norm * m_sigma / m_base_sigma; }
    Value intensity(Value b, Value l, Value T) {
        return std::exp(intensity_core(b, l, std::log(T)))
               * m_sigma / m_base_sigma;
    }
    ~IntensityInjector() {}
    Value next(Particle<Vector3, Value>& p) override {
        Value T = std::exp(this->m_rn.uniform_ab(std::log(m_Tmin),
                                                 std::log(m_Tmax)));
        Value b = this->m_rn.uniform_ab(m_bmin, m_bmax);
        Value l = this->m_rn.uniform_ab(m_lmin, m_lmax);
        init_Tbl(p, T, b, l, this->m_rn, this->m_t0, 1.1*this->m_radius);
        inject(p, this->m_radius);
        return intensity(b, l, T) * std::cos(b) * T * std::log(m_Tmax / m_Tmin)
               * (m_bmax - m_bmin) * (m_lmax - m_lmin) / getNorm();
    }
};

template<typename Vector3, typename Value>
class SampleInjector : public Injector<Vector3, Value>
{
private:
    hsize_t start;
    hsize_t end;
    hsize_t buf_idx;
    hsize_t buf_size;
    hsize_t mass_idx;
    std::string datafile;
    std::vector<Value> T;
    std::vector<Value> b;
    std::vector<Value> l;
    std::vector<Value> weight;
    hid_t fd;
    hid_t dxpl;
    hid_t fapl;
    hid_t ds_T;
    hid_t ds_bl;
    hid_t ds_w;
    void loadData() {
        hsize_t nrows = std::min(buf_size, end - start);

        if (nrows <= 0) { throw std::runtime_error("Run out of sample"); }

        T = read_dset<Value>(ds_T, {start}, {nrows}, dxpl);
        weight = read_dset<Value>(ds_w, {mass_idx, start}, {1, nrows}, dxpl);
        b = read_dset<Value>(ds_bl, {start, 0}, {nrows, 1}, dxpl);
        l = read_dset<Value>(ds_bl, {start, 1}, {nrows, 1}, dxpl);

        start += nrows;
        buf_idx = 0;
        buf_size = nrows;
    }
    void closeDsets() {
        H5Dclose(ds_T);
        H5Dclose(ds_w);
        H5Dclose(ds_bl);
    }
public:
    Value Tmin = std::numeric_limits<Value>::max();
    Value Tmax = 0;
    explicit SampleInjector(int seed=-1, Value t_radius=constants::rEarth)
        : buf_idx {0}, buf_size {10000} {
        this->m_radius = t_radius;
        this->setSeed(seed);
        dxpl = H5Pcreate(H5P_DATASET_XFER);
        fapl = H5Pcreate(H5P_FILE_ACCESS);
        datafile = "";
    }
    explicit SampleInjector(const std::string& filename,
                            unsigned flag=H5F_ACC_RDONLY,
                            int seed=-1, Value t_radius=constants::rEarth)
    : SampleInjector(seed, t_radius) {
        open(filename, flag);
    }
    SampleInjector(const SampleInjector&) = delete;
    SampleInjector& operator=(const SampleInjector&) = delete;
    ~SampleInjector() {
        H5Pclose(dxpl);
        H5Pclose(fapl);
        if (!datafile.empty()) {
            closeDsets();
            H5Fclose(fd);
        }
    }
    void open(const std::string& filename, unsigned int flag=H5F_ACC_RDONLY) {
        if (datafile == filename) { return; }
        if (!datafile.empty()) { closeDsets(); H5Fclose(fd); }
        if (H5Fis_hdf5(filename.c_str())) {
            datafile = filename;
            fd = H5Fopen(filename.c_str(), flag, fapl);
            ds_T = H5Dopen(fd, "ran_E", H5P_DEFAULT);
            ds_w = H5Dopen(fd, "weight", H5P_DEFAULT);
            ds_bl = H5Dopen(fd, "ran_bl", H5P_DEFAULT);
        } else {
            throw std::invalid_argument(filename + " is not a hdf5 file");
        }
    }
    void setBufferSize(hsize_t sz) { buf_size = sz; }
    void setDataRange(hsize_t t_start, hsize_t t_end) {
        start = t_start;
        end = t_end;
    }
    void setWeightMassIdx(hsize_t idx) {
        mass_idx = idx;
    }
    Value next(Particle<Vector3, Value>& p) override {
        if (buf_idx == buf_size || T.size() == 0) {
            loadData();
        }
        Value w = weight[buf_idx];
        Value Ti = T[buf_idx];
        Tmin = std::min(Ti, Tmin);
        Tmax = std::max(Ti, Tmax);
        // TODO other arguments
        init_Tbl(p, Ti, b[buf_idx], l[buf_idx], this->m_rn);
        ++buf_idx;
        inject(p, this->m_radius);
        return w;
    }
};


template<typename Vector3=Vector3d<double>, typename Value=double>
class SourceInjector : public Injector<Vector3, Value>
{
private:
    std::string method;
    GSLInterpolator<Value> inverse_cdf_r;
    GSLInterpolator<Value> inverse_cdf_T;
public:
    explicit SourceInjector(int t_seed=-1) { this->setSeed(t_seed); }
    SourceInjector(const std::vector<Value>& t_xi_T,
                   const std::vector<Value>& t_inv_cdf_T,
                   const std::vector<Value>& t_xi_r,
                   const std::vector<Value>& t_inv_cdf_r,
                   const std::string& t_method="linear",
                   int t_seed=-1) {
        this->method = t_method;
        this->setSeed(t_seed);
        interpolate_r(t_xi_r, t_inv_cdf_r);
        interpolate_T(t_xi_T, t_inv_cdf_T);
    }
    SourceInjector(const std::string& t_r_file,
                   const std::string& t_T_file,
                   const std::string& t_method="linear",
                   int seed=-1) {
        this->method = t_method;
        this->setSeed(seed);
        std::vector<Value> xi_r;
        std::vector<Value> xi_T;
        std::vector<Value> inv_cdf_r;
        std::vector<Value> inv_cdf_T;
        std::ifstream ifs;

        ifs.open(t_r_file);
        Value tmp;
        while (ifs >> tmp) {
            xi_r.push_back(tmp);
            ifs >> tmp;
            inv_cdf_r.push_back(tmp * units::cm);
        }
        ifs.close();

        ifs.open(t_T_file);
        while (ifs >> tmp) {
            xi_T.push_back(tmp);
            ifs >> tmp;
            inv_cdf_T.push_back(tmp * units::MeV);
        }
        ifs.close();

        interpolate_r(xi_r, inv_cdf_r);
        interpolate_T(xi_T, inv_cdf_T);
    }
    void interpolate_r(const std::vector<Value>& t_r,
                       const std::vector<Value>& t_spectrum) {
        inverse_cdf_r.interpolate(t_r, t_spectrum, this->method);
    }
    void interpolate_T(const std::vector<Value>& t_T,
                       const std::vector<Value>& t_spectrum) {
        inverse_cdf_T.interpolate(t_T, t_spectrum, this->method);
    }
    ~SourceInjector() {}
    Value next(Particle<Vector3, Value>& p) override {
        Value T = inverse_cdf_T(this->m_rn.uniform_xi());
        Value r = inverse_cdf_r(this->m_rn.uniform_xi());
        // TODO check the boundary
        init_rT_isotropic(p, r, T, this->m_rn, this->m_t0);
        return 1.0;
    }
};

} // namespace darkprop
