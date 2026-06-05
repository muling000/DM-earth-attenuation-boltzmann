/**
 * @file HDF5 IO.
 * Parallel hdf5 is supported.
 */
#pragma once

#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <type_traits>
#include <fstream>
#include <hdf5.h>
#include "Base.hpp"

namespace darkprop {

/**
 * Determine HDF5 hid_t type.
 */
template<typename T>
constexpr hid_t hid_type(T var=T()) {
    if constexpr (std::is_same<T, int>::value) {
        return H5T_STD_I32LE;
    }
    if constexpr (std::is_same<T, long>::value) {
        return H5T_STD_I64LE;
    }
    if constexpr (std::is_same<T, unsigned int>::value) {
        return H5T_STD_U32LE;
    }
    if constexpr (std::is_same<T, unsigned long>::value) {
        return H5T_STD_U64LE;
    }
    if constexpr (std::is_same<T, const char *>::value) {
        hid_t type = H5Tcopy(H5T_C_S1);
        if (var) {
            H5Tset_size(type, strlen(var));
        }
        return type;
    }
    if constexpr (std::is_same<T, float>::value) {
        return H5T_IEEE_F32LE;
    }
    return H5T_IEEE_F64LE;
}

template<typename T>
void write_attr(hid_t& dset, const char *attr_name, T attr_value)
{
    if (attr_name) {
        hid_t aspace = H5Screate_simple(0, NULL, NULL);
        hid_t type = hid_type(attr_value);
        hid_t attr = H5Acreate(dset, attr_name, type, aspace,
                            H5P_DEFAULT, H5P_DEFAULT);
        if constexpr (std::is_same<T, const char *>::value) {
            H5Awrite(attr, type, attr_value);
            H5Tclose(type);
        } else {
            H5Awrite(attr, type, &attr_value);
        }
        H5Aclose(attr);
        H5Sclose(aspace);
    }
}

template<typename T>
void write_attr(hid_t& group, const char *dset_name, const char *attr_name, T attr_value)
{
    hid_t dset = H5Dopen(group, dset_name, H5P_DEFAULT);
    write_attr(dset, attr_name, attr_value);
    H5Dclose(dset);
}

template<typename D, typename A>
void write_scalar(hid_t& group,
                  const char *dset_name,
                  D dset_value,
                  const char *attr_name,
                  A attr_value)
{
    hid_t type = hid_type(dset_value);
    hid_t dspace = H5Screate_simple(0, NULL, NULL);
    hid_t dset = H5Dcreate(group, dset_name, type, dspace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dset_value);

    if (attr_name) {
        write_attr(dset, attr_name, attr_value);
    }

    H5Dclose(dset);
    H5Sclose(dspace);
}

template<typename D, typename A>
void write_scalar(const char *filename,
                  const char *dset_name,
                  D dset_value,
                  const char *attr_name,
                  A attr_value)
{
    hid_t file;
    if (std::filesystem::exists(filename)) {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    write_scalar(file, dset_name, dset_value, attr_name, attr_value);
    H5Fclose(file);
}

template<typename D, typename A>
void write_scalar(const char *filename,
                  const char *groupname,
                  const char *dset_name,
                  D dset_value,
                  const char *attr_name,
                  A attr_value)
{
    hid_t file, group;
    if (std::filesystem::exists(filename)) {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Lexists(file, groupname, H5P_DEFAULT)) {
        group = H5Gopen(file, groupname, H5P_DEFAULT);
    } else {
        group = H5Gcreate(file, groupname,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    std::string data_path = std::string(groupname) + "/" + dset_name;
    if (!H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        write_scalar(group, dset_name, dset_value, attr_name, attr_value);
    }
    H5Gclose(group);
    H5Fclose(file);
}

template<typename D>
void write_scalar(const char *filename, const char *dset_name, D value)
{
    write_scalar(filename, dset_name, value, nullptr, nullptr);
}

/**
 * Extend at the 0th dimension
 */
template<typename C>
void append_array(hid_t& group, const char *dset_name, const C& new_data)
{
    hid_t dset = H5Dopen(group, dset_name, H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    std::vector<hsize_t> oldims(H5Sget_simple_extent_ndims(dspace));
    int ndim = H5Sget_simple_extent_dims(dspace, oldims.data(), NULL);
    std::vector<hsize_t> dims(oldims);
    dims[0] = new_data.size();
    for (std::size_t i = 1; i < oldims.size(); ++i) {
        dims[0] /= oldims[i];
    }
    std::vector<hsize_t> extdims(oldims);
    H5Sclose(dspace);
    extdims[0] = oldims[0] + dims[0];
    H5Dset_extent(dset, extdims.data());
    dspace = H5Dget_space(dset);
    std::vector<hsize_t> start(ndim, 0);
    start[0] = oldims[0];
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET,
                        start.data(), NULL, dims.data(), NULL);
    hid_t memspace = H5Screate_simple(ndim, dims.data(), NULL);
    H5Dwrite(dset, H5Dget_type(dset), memspace, dspace,
             H5P_DEFAULT, new_data.data());
    H5Sclose(memspace);
    H5Sclose(dspace);
    H5Dclose(dset);
}

template<typename C>
void write_array(hid_t& group,
                 const char *dset_name,
                 const C& data,
                 const std::vector<hsize_t>& dims,
                 hid_t& type,
                 hsize_t chunk_size)
{
    std::vector<hsize_t> maxdims(dims);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (chunk_size > 0) {
        maxdims[0] = H5S_UNLIMITED;
        std::vector<hsize_t> chunk(dims);
        chunk[0] = chunk_size;
        H5Pset_chunk(dcpl, dims.size(), chunk.data());
        H5Pset_deflate(dcpl, 4);
    }

    hid_t dspace = H5Screate_simple(dims.size(), dims.data(), maxdims.data());
    hid_t dset = H5Dcreate(group, dset_name, type, dspace,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(dset);
    H5Sclose(dspace);
    H5Pclose(dcpl);
}


/**
 * Dump events into a hdf5 file using galactic coordinates b and l.
 * This function is used to store large data. The stored data type is the same
 * as \p R and \p m. Please refer to hdfgroup to determine the chunk size, the
 * default is 10000.
 */
template<typename TStore=float, typename Vector3, typename Value>
void h5dump_events_galactic(const char* filename,
                            const char* groupname,
                            const std::vector<Event<Vector3, Value>>& events,
                            Value R,  ///< radius [km]
                            Value m,  ///< mass [GeV]
                            hsize_t chunk_size=10000)
{
    hid_t file, group;
    std::string data_path;

    if (std::filesystem::exists(filename)) {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (H5Lexists(file, groupname, H5P_DEFAULT)) {
        group = H5Gopen(file, groupname, H5P_DEFAULT);
    } else {
        group = H5Gcreate(file, groupname,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    // dark matter mass
    data_path = std::string(groupname) + "/mchi";
    if (!H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        write_scalar(group, "mchi", m, "unit", "GeV");
    }

    // Radius
    data_path = std::string(groupname) + "/R";
    if (!H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        write_scalar(group, "R", R, "unit", "km");
    }

    std::size_t N = events.size();
    std::vector<TStore> Ts(N);
    std::vector<TStore> rs(N*2);
    std::vector<TStore> ps(N*2);
    for (std::size_t i = 0; i != N; ++i) {
        const auto& e = events[i];
        Ts[i] = static_cast<TStore>(e.T);
        rs[2*i] = static_cast<TStore>(std::asin(e.r[2] / e.r.norm()));
        rs[2*i+1] = static_cast<TStore>(std::atan2(e.r[1], e.r[0]));
        ps[2*i] = static_cast<TStore>(std::asin(e.p3[2] / e.p3.norm()));
        ps[2*i+1] = static_cast<TStore>(std::atan2(e.p3[1], e.p3[0]));
    }

    hid_t type = hid_type<TStore>();

    // kinetic energy
    data_path = std::string(groupname) + "/T";
    if (H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        append_array(group, "T", Ts);
    } else {
        write_array(group, "T", Ts, {N}, type, chunk_size);
        write_attr(group, "T", "unit", "GeV");
    }

    // normalized position r (b, l)
    data_path = std::string(groupname) + "/r";
    if (H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        append_array(group, "r", rs);
    } else {
        write_array(group, "r", rs, {N, 2}, type, chunk_size);
        write_attr(group, "r", "unit", "rad");
        write_attr(group, "r", "dimension", "(b, l)");
    }

    // normalized 3 momentum p3 (b, l)
    data_path = std::string(groupname) + "/p";
    if (H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        append_array(group, "p", ps);
    } else {
        write_array(group, "p", ps, {N, 2}, type, chunk_size);
        write_attr(group, "p", "unit", "rad");
        write_attr(group, "p", "dimension", "(b, l)");
    }

    H5Gclose(group);
    H5Fclose(file);
}

template<typename TStore=float, typename Vector3, typename Value>
void h5dump_events_galactic(const std::string& filename,
                            const std::string& groupname,
                            const std::vector<Event<Vector3, Value>>& events,
                            Value R,  ///< radius [km]
                            Value m,  ///< mass [GeV]
                            hsize_t chunk_size=10000)
{
    h5dump_events_galactic<TStore>(filename.c_str(), groupname.c_str(),
                                   events, R, m, chunk_size);
}


template<typename arr>
void h5dump_array(const char *filename,
                  const char *groupname,
                  const char *dset_name,
                  const arr& v,
                  const std::vector<hsize_t>& dims,
                  hsize_t chunk_size)
{
    hid_t file, group, type = hid_type(v[0]);

    if (std::filesystem::exists(filename)) {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (H5Lexists(file, groupname, H5P_DEFAULT)) {
        group = H5Gopen(file, groupname, H5P_DEFAULT);
    } else {
        group = H5Gcreate(file, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    std::string data_path = std::string(groupname) + "/" + dset_name;
    if (H5Lexists(file, data_path.c_str(), H5P_DEFAULT)) {
        append_array(group, dset_name, v);
    } else {
        write_array(group, dset_name, v, dims, type, chunk_size);
    }

    H5Gclose(group);
    H5Fclose(file);
}

template<typename arr>
void h5dump_array(const std::string& filename,
                  const std::string& groupname,
                  const std::string& dset_name,
                  const arr& v,
                  const std::vector<hsize_t>& dims,
                  hsize_t chunk_size)
{
    h5dump_array(filename.c_str(), groupname.c_str(), dset_name.c_str(),
                 v, dims, chunk_size);
}

/**
 * Dump events into a hdf5 file using Cartesian coordinates.
 * This function is used to dump small data at once. Complete but also
 * redundant information of events is stored.
 * The stored data type is \p Value. Please refer to hdfgroup to determine the
 * chunk size, the default is not chunked (\p chunk_size = 0).
 */
template<typename TStore=float, typename Vector3, typename Value>
void h5dump_events(const char* filename,
                   const char* groupname,
                   const std::vector<Event<Vector3, Value>>& events,
                   hsize_t chunk_size=0)
{
    hid_t file, group, type = hid_type<TStore>();

    if (std::filesystem::exists(filename)) {
        file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Lexists(file, groupname, H5P_DEFAULT)) {
        group = H5Gopen(file, groupname, H5P_DEFAULT);
    } else {
        group = H5Gcreate(file, groupname,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    size_t N = events.size();
    std::vector<TStore> ts(N);
    std::vector<TStore> Ts(N);
    std::vector<TStore> rs(3*N);
    std::vector<TStore> ps(3*N);
    for (size_t i = 0; i != N; ++i) {
        const auto& e = events[i];
        ts[i] = e.t;
        Ts[i] = e.T;
        rs[3*i + 0] = e.r[0];
        rs[3*i + 1] = e.r[1];
        rs[3*i + 2] = e.r[2];
        ps[3*i + 0] = e.p3[0];
        ps[3*i + 1] = e.p3[1];
        ps[3*i + 2] = e.p3[2];
    }

    write_array(group, "t", ts, {N}, type, chunk_size);
    write_attr(group, "t", "unit", "s");
    write_array(group, "T", Ts, {N}, type, chunk_size);
    write_attr(group, "T", "unit", "GeV");
    write_array(group, "r", rs, {N, 3}, type, chunk_size);
    write_attr(group, "r", "unit", "km");
    write_array(group, "p", ps, {N, 3}, type, chunk_size);
    write_attr(group, "p", "unit", "GeV");
}


template<typename Value>
void h5load_array(const char *filename,
                  const char *groupname,
                  const char *dset_name,
                  std::vector<Value>& v,
                  std::vector<hsize_t>& dims,
                  const hid_t& fapl=H5P_DEFAULT)
{
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, fapl);
    hid_t group = H5Gopen(file, groupname, H5P_DEFAULT);
    hid_t dset = H5Dopen(group, dset_name, H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    hid_t type = hid_type<Value>();

    int ndims = H5Sget_simple_extent_ndims(dspace);
    dims.resize(ndims);
    H5Sget_simple_extent_dims(dspace, dims.data(), NULL);

    std::size_t data_size = dims[0];
    for (std::size_t i = 1; i < dims.size(); ++i) {
        data_size *= dims[i];
    }
    v.resize(data_size);

    H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Gclose(group);
    H5Fclose(file);
}

template<typename Value>
Value h5load_scalar(const char *filename,
                    const char *groupname,
                    const char *dset_name,
                    const hid_t& fapl=H5P_DEFAULT)
{
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, fapl);
    hid_t group = H5Gopen(file, groupname, H5P_DEFAULT);
    hid_t dset = H5Dopen(group, dset_name, H5P_DEFAULT);
    hid_t type = hid_type<Value>();

    Value value;
    H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

    H5Dclose(dset);
    H5Gclose(group);
    H5Fclose(file);
    return value;
}

template<typename Value>
void merge_dataset(const char *tofilename,
                   const char *fromfilename,
                   const char *groupname,
                   const char *dset_name,
                   hsize_t chunk_size=10000)
{
    std::vector<Value> v;
    std::vector<hsize_t> dims;
    h5load_array(fromfilename, groupname, dset_name, v, dims);
    h5dump_array(tofilename, groupname, dset_name, v, dims, chunk_size);

}

template<typename Value>
void merge_dataset(const std::string& tofilename,
                   const std::string& fromfilename,
                   const std::string& groupname,
                   const std::string& dset_name)
{
    merge_dataset<Value>(tofilename.c_str(), fromfilename.c_str(),
                         groupname.c_str(), dset_name.c_str());
}


inline hid_t make_group(hid_t loc_id, const std::string& gname) {
    return H5Gcreate(loc_id, gname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}


template<typename TStore=float>
hid_t make_dset(hid_t loc_id, const std::string& dsetname,
                std::vector<hsize_t> cur_dims={},
                std::vector<hsize_t> max_dims={},
                std::vector<hsize_t> chunk={})
{
    hid_t type = hid_type<TStore>();
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    TStore fill_value = 0;
    if constexpr (std::is_floating_point<TStore>::value) {
        fill_value = NAN;
    }
    H5Pset_fill_value(dcpl, type, &fill_value);
    if (chunk.size() != 0) {
        H5Pset_chunk(dcpl, chunk.size(), chunk.data());
    }
    hid_t dspace = H5Screate_simple(static_cast<int>(cur_dims.size()),
                                    cur_dims.data(), max_dims.data());
    hid_t ds = H5Dcreate(loc_id, dsetname.c_str(), type, dspace,
                         H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Sclose(dspace);
    H5Pclose(dcpl);
    return ds;
}

template<typename TStore=float, typename C>
void write_dset(hid_t ds, const C& data,
                std::vector<hsize_t> offset,
                std::vector<hsize_t> count,
                hid_t xferpl=H5P_DEFAULT)
{
    hid_t memsp = H5Screate_simple(count.size(), count.data(), NULL);
    hid_t filesp = H5Dget_space(ds);
    H5Sselect_hyperslab(filesp, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);
    H5Dwrite(ds, hid_type<TStore>(), memsp, filesp, xferpl, data.data());
    H5Sclose(memsp);
    H5Sclose(filesp);
}

template<typename TStore=float>
inline void write_dset(hid_t ds, TStore data, hid_t xferpl=H5P_DEFAULT)
{
    // TODO check size
    H5Dwrite(ds, hid_type<TStore>(), H5S_ALL, H5S_ALL, xferpl, &data);
}

template<typename TStore=float, typename Value>
std::map<std::string, hid_t> make_dsets(
    hid_t fid,
    const std::vector<std::string>& depth_groups,
    const std::vector<std::string>& ring_groups,
    const std::vector<Value>& depth,
    const std::vector<Value>& ring_bins,
    Value radius,
    std::size_t sample_sz,
    std::size_t chunk_size,
    bool store_time=false)
{
    std::map<std::string, hid_t> dsets;
    std::string dsetname = "depth";
    // TODO automatic conversion
    std::vector<TStore> depth_store(depth.size());
    for (std::size_t di = 0; di != depth.size(); ++di) {
        depth_store[di] = static_cast<TStore>(depth[di] / units::km);
    }
    dsets[dsetname] = make_dset<TStore>(fid, dsetname, {depth_store.size()});
    // TODO collective
    write_dset<TStore>(dsets[dsetname], depth_store, {0}, {depth.size()});
    write_attr(dsets[dsetname], "unit", "km");
    write_attr(dsets[dsetname], "description",
               "This array corresponds to the depth_{i} groups.");
    for (std::size_t di = 0; di != depth_groups.size(); ++di) {
        hid_t gid = make_group(fid, depth_groups[di]);
        H5Gclose(gid);
        dsetname = depth_groups[di] + "/depth";
        dsets[dsetname] = make_dset<TStore>(fid, dsetname);
        write_dset<TStore>(dsets[dsetname], depth[di] / units::km);
        write_attr(dsets[dsetname], "unit", "km");

        dsetname = depth_groups[di] + "/R";
        dsets[dsetname] = make_dset<TStore>(fid, dsetname);
        write_dset<TStore>(dsets[dsetname], (radius - depth[di]) / units::km);
        write_attr(dsets[dsetname], "unit", "km");

        if (ring_bins.size() > 2) {
            dsetname = depth_groups[di] + "/ring_bins";
            std::vector<TStore> bins_store(ring_bins.begin(), ring_bins.end());
            dsets[dsetname] = make_dset<TStore>(fid, dsetname,
                                                {bins_store.size()});
            write_dset(dsets[dsetname], bins_store, {0}, {bins_store.size()});
            write_attr(dsets[dsetname], "unit", "rad");

        }

        for (const auto& ring_group : ring_groups) {
            std::string group = depth_groups[di] + "/";
            if (ring_groups.size() > 1) {
                group = group + ring_group + "/";
                gid = make_group(fid, group);
                H5Gclose(gid);
            }
            dsetname = group + "T";
            dsets[dsetname] = make_dset<TStore>(fid, dsetname, {sample_sz},
                                                {H5S_UNLIMITED}, {chunk_size});
            write_attr(dsets[dsetname], "unit", "GeV");

            if (store_time) {
                dsetname = group + "t";
                dsets[dsetname] = make_dset<TStore>(fid, dsetname, {sample_sz},
                                                    {H5S_UNLIMITED},
                                                    {chunk_size});
                write_attr(dsets[dsetname], "unit", "second");
            }

            dsetname = group + "weight";
            dsets[dsetname] = make_dset<TStore>(fid, dsetname, {sample_sz},
                                                {H5S_UNLIMITED}, {chunk_size});

            dsetname = group + "r";
            dsets[dsetname] = make_dset<TStore>(fid, dsetname, {sample_sz, 2},
                                        {H5S_UNLIMITED, 2}, {chunk_size, 2});
            write_attr(dsets[dsetname], "unit", "rad");

            dsetname = group + "p";
            dsets[dsetname] = make_dset<TStore>(fid, dsetname, {sample_sz, 2},
                                        {H5S_UNLIMITED, 2}, {chunk_size, 2});
            write_attr(dsets[dsetname], "unit", "rad");

            dsetname = group + "Nsim";
            dsets[dsetname] = make_dset<std::size_t>(fid, dsetname);
            dsetname = group + "Nsample";
            dsets[dsetname] = make_dset<std::size_t>(fid, dsetname);
        }
    }
    return dsets;
}

inline void close_dsets(std::map<std::string, hid_t>& dsets)
{
    for (auto it = dsets.begin(); it != dsets.end(); ++it) {
        H5Dclose(it->second);
    }
}

template<typename TStore=float, typename Vector3, typename Value>
void h5dump_events_galactic(std::map<std::string, hid_t> dsets,
                            const std::string& groupname,
                            const std::vector<Event<Vector3, Value>>& events,
                            std::size_t offset,
                            bool store_time=false)
{
    std::size_t N = events.size();
    std::vector<TStore> Ts(N);
    std::vector<TStore> rs(N*2);
    std::vector<TStore> ps(N*2);
    for (std::size_t i = 0; i != N; ++i) {
        Ts[i] = static_cast<TStore>(events[i].T);
        rs[2*i] = static_cast<TStore>(std::asin(events[i].r[2]
                                                / events[i].r.norm()));
        Value l = std::atan2(events[i].r[1], events[i].r[0]);
        rs[2*i+1] = static_cast<TStore>(l >= 0 ? l : l + 2.0 * M_PI);
        ps[2*i] = static_cast<TStore>(std::asin(events[i].p3[2]
                                                / events[i].p3.norm()));
        l = std::atan2(events[i].p3[1], events[i].p3[0]);
        ps[2*i+1] = static_cast<TStore>(l >= 0 ? l : l + 2.0 * M_PI);
    }
    write_dset<TStore>(dsets[groupname + "/T"], Ts, {offset}, {N});
    write_dset<TStore>(dsets[groupname + "/r"], rs, {offset, 0}, {N, 2});
    write_dset<TStore>(dsets[groupname + "/p"], ps, {offset, 0}, {N, 2});
    if (store_time) {
        std::vector<TStore> ts(N);
        for (std::size_t i = 0; i != N; ++i) {
            ts[i] = static_cast<TStore>(events[i].t);
        }
        write_dset(dsets[groupname + "/t"], ts, {offset}, {N});
    }
}

template<typename T>
std::vector<T> read_dset(hid_t ds,
                         const std::vector<hsize_t>& offset={},
                         const std::vector<hsize_t>& count={},
                         hid_t xferpl=H5P_DEFAULT)
{
    hid_t type = hid_type<T>();
    hid_t filespace = H5Dget_space(ds);
    if (count.size() != 0) {
        hsize_t N = 1;
        for (auto i : count) { N *= i; }
        std::vector<T> data(N);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), NULL,
                            count.data(), NULL);
        hid_t memspace = H5Screate_simple(1, &N, NULL);
        H5Dread(ds, type, memspace, filespace, xferpl, data.data());
        H5Sclose(memspace);
        H5Sclose(filespace);
        return data;
    }
    int ndim = H5Sget_simple_extent_ndims(filespace);
    std::vector<hsize_t> dims(ndim);
    H5Sget_simple_extent_dims(filespace, dims.data(), NULL);
    hsize_t N = 1;
    for (auto i : dims) { N *= i; }
    std::vector<T> data(N);
    hid_t memspace = H5Screate_simple(1, &N, NULL);
    H5Dread(ds, type, memspace, H5S_ALL, xferpl, data.data());
    H5Sclose(memspace);
    H5Sclose(filespace);
    return data;
}

template<typename T>
std::vector<T> read_dset(hid_t fd, const std::string& dsetname,
                         const std::vector<hsize_t>& offset={},
                         const std::vector<hsize_t>& count={},
                         hid_t xferpl=H5P_DEFAULT)
{
    // TODO check property list dapl
    hid_t ds = H5Dopen(fd, dsetname.c_str(), H5P_DEFAULT);
    auto data = read_dset<T>(ds, offset, count, xferpl);
    H5Dclose(ds);
    return data;
}

/**
 * Read dataset from a hdf5 file.
 * Return a one dimensional \p std::vector of type \p T flattened from the data read.
 *
 * The \p offset and \p count parameters are used to specify the part of the data set to be
 * read, if they are empty, read the entire dataset.
 */
template<typename T>
std::vector<T> read_dset(
    const std::string& filename,            ///< hdf5 filename
    const std::string& dsetname,            ///< dataset name (full path)
    const std::vector<hsize_t>& offset={},  ///< offset
    const std::vector<hsize_t>& count={},   ///< count
    unsigned int flag=H5F_ACC_RDONLY,       ///< flag pass to H5Fopen
    hid_t fapl=H5P_DEFAULT,                 ///< fapl pass to H5Fopen
    hid_t xferpl=H5P_DEFAULT)               ///< xferpl pass to H5Dread
{
    hid_t fd = H5Fopen(filename.c_str(), flag, fapl);
    auto data = read_dset<T>(fd, dsetname, offset, count, xferpl);
    H5Fclose(fd);
    return data;
}

template<typename Vector3, typename Value>
void dump_track_csv(const std::vector<Event<Vector3, Value>>& events,
                    const std::string& filename)
{
    std::ofstream ofs(filename);
    ofs << "# r_x [km] r_y [km] r_z [km] p_x [GeV] p_y [GeV] p_z [GeV] T [GeV] t [sec]\n";
    ofs.setf(std::ios_base::scientific);
    ofs.precision(16);
    for (const auto& event : events) {
        ofs << event.r[0] << " " << event.r[1] << " " << event.r[2] << " "
            << event.p3[0] << " " << event.p3[1] << " " << event.p3[2] << " "
            << event.T << " " << event.t << '\n';
    }
}

} // namespace darkprop
