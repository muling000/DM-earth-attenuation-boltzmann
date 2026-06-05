#pragma once
#include <fstream>
#include <initializer_list>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <valarray>
#include <variant>
#include <vector>

namespace yuc {
struct config;
namespace __config_detail {
using nil_t = std::monostate;
using arr_t = std::vector<config>;
using obj_t = std::vector<std::pair<std::string, config>>;
using str_t = std::string;
using dbl_t = double;
using int_t = long;
using bol_t = bool;

struct time_type {
public:
  short year = 0, month = 0, week = 0, day = 0;
  short hour = 0, min = 0, offset_hour = 0, offset_min = 0;
  double sec = 0;
  friend class config;

public:
  [[deprecated("detailed datetime not implemented")]] double
  to_utc(void) const {
    return 0;
  }
  template <typename T> bool operator==(const T &) const { return false; }
};

using tim_t = time_type;

using var_t =
    std::variant<nil_t, arr_t, obj_t, str_t, dbl_t, int_t, bol_t, tim_t>;
inline constexpr size_t var_n = std::variant_size_v<var_t>;

template <typename T, size_t I = 0>
inline constexpr size_t var_i =
    std::is_same_v<std::remove_cv_t<T>, std::variant_alternative_t<I, var_t>>
        ? I
        : var_i<T, I + 1>;
template <typename T> inline constexpr size_t var_i<T, var_n> = var_n;

template <typename> struct block_implicit : std::false_type {};
template <typename T>
struct block_implicit<std::initializer_list<T>> : std::true_type {};
template <typename T> struct block_implicit<std::vector<T>> : std::true_type {};
template <typename T>
struct block_implicit<std::valarray<T>> : std::true_type {};
template <> struct block_implicit<long> : std::true_type {};

template <typename T, size_t I = 0>
inline constexpr bool from_union =
    std::is_same_v<std::remove_cv_t<T>, std::variant_alternative_t<I, var_t>>
        ? true
        : from_union<T, I + 1>;
template <typename T> inline constexpr bool from_union<T, var_n> = false;

template <typename T, size_t I = 0, size_t Candidate = 0>
inline T config_convert_rt(const var_t &v) {
  static_assert(I < var_n || Candidate, "conversion is impossible");
  if constexpr (I == var_n) {
    throw std::logic_error("unable to convert: " +
                           std::string(__PRETTY_FUNCTION__));
  } else if constexpr (std::is_convertible_v<
                           std::variant_alternative_t<I, var_t>, T>) {
    if (v.index() == I) {
      return std::get<I>(v);
    }
    return config_convert_rt<T, I + 1, Candidate + 1>(v);
  } else {
    return config_convert_rt<T, I + 1, Candidate>(v);
  }
}

inline config &get_context(config &root, std::istream &is, char del,
                           const std::string &fn);

// ---------- string utilities ----------

inline bool peek_forward(std::istream &is, const std::string_view &expect) {
  auto pos = is.tellg();
  for (auto c : expect) {
    if (c != is.get()) {
      return is.clear(), is.seekg(pos), false;
    }
  }
  return is.clear(), is.seekg(pos), true;
}

inline std::string string_unescape(const std::string_view &ess) {
  auto iti = ess.cbegin(), ite = ess.cend();
  std::string raw(ess);
  while (iti != ite && *iti++ != '\\')
    ;
  if (iti == ite) { // no escape characters in most calls
    return raw;
  }
  auto ito = raw.begin() + (--iti - ess.begin());
  while (iti != ite) {
    if ((*ito = *iti++) == '\\' && iti != ite) {
      *iti == '\\'  ? *ito = '\\'
      : *iti == '"' ? *ito = '"'
      : *iti == 'a' ? *ito = '\a'
      : *iti == 'b' ? *ito = '\b'
      : *iti == 'f' ? *ito = '\f'
      : *iti == 'n' ? *ito = '\n'
      : *iti == 'r' ? *ito = '\r'
      : *iti == 't' ? *ito = '\t'
      : *iti == 'v' ? *ito = '\v'
                    : *(++ito) = *iti;
      ++iti;
    }
    ++ito;
  }
  return std::string(raw.begin(), ito);
}

inline std::string string_escape(const std::string_view &sv) {
  std::string res;
  res.reserve(sv.size());
  for (const char &c : sv) {
    res += c == '\\'         ? "\\\\"
           : c == '"'        ? "\\\""
           : c == '\a'       ? "\\a"
           : c == '\b'       ? "\\b"
           : c == '\f'       ? "\\f"
           : c == '\n'       ? "\\n"
           : c == '\r'       ? "\\r"
           : c == '\t'       ? "\\t"
           : c == '\v'       ? "\\v"
           : std::isprint(c) ? (res.push_back(c), "")
                             : dynamic_cast<std::ostringstream &&>(
                                   std::ostringstream()
                                   << (uint8_t(c) > '\xf' ? "\\x" : "\\x0")
                                   << std::hex << (int)c)
                                   .str()
                                   .data();
  }
  return res;
}

// escaped @quote_char will NOT be unescaped in @dest
// next unescaped @quote_char is discard but not read into @dest
// caller set @quote_end_char and @escape_char the same to not escape anything
// failbit is set if stream ends prematurely
inline auto read_quoted(std::istream &is, char end_quote_c = '"',
                        char escape_c = '\\') {
  std::string temp;
  char c;
  while (is.get(c)) {
    if (c == end_quote_c) {
      break;
    } else if (c == escape_c) {
      temp.push_back(c), is.get(c);
    }
    temp.push_back(c);
  }
  return temp;
}
inline auto &read_quoted(std::istream &is, std::string &dest,
                         char end_quote_c = '"', char escape_c = '\\') {
  auto temp = read_quoted(is, end_quote_c, escape_c);
  return is.good() ? dest = temp, is : is;
}

inline std::ostream &write_quoted(std::ostream &os,
                                  const std::string_view &sv) {
  return os << '"' << string_escape(sv) << '"';
}
inline std::string string_quote(const std::string_view &sv) {
  std::string res = "\"";
  res += string_escape(sv);
  res += '"';
  return res;
}

// start from unquoted pos, stop at unquoted space or delimiter if given
// any unquoted/unescaped quote mark is removed
// eg: ab"c d' ef\"g\\"hi'jk"lm\'n
// --> abc d' ef"g\hijk"lm\n
inline auto read_word(std::istream &is, char del = EOF) {
  std::string temp;
  is >> std::ws;
  while (true) {
    if (char c; is.get(c).fail()) {
      if (temp.size()) { // clear failbit read something
        is.clear();
      }
      return temp;
    } else if (::isspace(c) || c == del) {
      is.unget();
      return temp;
    } else if (c == '\'') {
      temp += read_quoted(is, c, c);
    } else if (c == '"') {
      temp += string_unescape(read_quoted(is));
    } else if (c == '\'') {
      temp += read_quoted(is, '\'', EOF);
    } else {
      temp.push_back(c);
    }
  }
}
inline auto &read_word(std::istream &is, std::string &dest, char del = EOF) {
  auto temp = read_word(is, del);
  return is.good() ? dest = temp, is : is;
}

inline auto string_split(const std::string_view &sv,
                         const std::string_view &sep = "",
                         bool skip_empty = true) {
  std::vector<std::string_view> vec;
  const auto nsp = sep.size();
  const auto ite = sv.end();
  auto it0 = sv.begin(), it1 = it0;
  while (it1 < ite) {
    if (nsp == 0) { // use white spaces as seperater
      it1 = std::find_if(it0, ite, ::isspace);
    } else {
      for (it1 = it0;
           it1 < ite && (it1 > ite - nsp || std::string_view(it1, nsp) != sep);
           ++it1) {
      }
    }
    if (it1 > it0 || !skip_empty) {
      vec.emplace_back(it0, it1 - it0);
    }
    it0 = it1 + nsp + !nsp;
  }
  return vec;
}

inline std::istream &stream_trim(std::istream &is,
                                 const std::string_view &cmt = "") {
  while (is >> std::ws, cmt.size() && is.peek() == cmt[0] && is.good()) {
    const auto pos = is.tellg();
    for (const char c : cmt) {
      if (c != is.get()) {
        return is.clear(), is.seekg(pos), is;
      }
    }
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return is;
}

inline std::string_view string_trim(const std::string_view &sv, int dir = 0) {
  auto start =
      dir < 0 ? sv.begin() : std::find_if_not(sv.begin(), sv.end(), ::isspace);
  auto stop =
      dir > 0 || start == sv.end()
          ? sv.end()
          : std::find_if_not(sv.rbegin(), std::make_reverse_iterator(start),
                             ::isspace)
                .base();
  return {start, static_cast<size_t>(stop - start)};
}

// ----------- io exception ------------

class parse_error : public std::runtime_error {
public:
  std::string fn, msg, line;
  size_t ln, cn;

  explicit parse_error(std::istream &is, const std::string &filename,
                       const std::string &message)
      : runtime_error(""), fn(filename), msg(message), ln(1), cn(1) {
    bool eof = is.eof();
    const auto bits = is.rdstate();
    is.clear();
    const auto cpos = is.tellg();
    auto lpos = cpos;
    do {
      while (is.unget() && is.peek() != '\n') {
        ++cn, lpos = is.tellg();
      }
    } while (cn == 1 && eof && (eof = false, true));
    do {
      ln += is.peek() == '\n';
    } while (is.unget());
    is.clear();
    is.seekg(lpos);
    std::getline(is.seekg(lpos), line);
    is.seekg(cpos), is.setstate(bits); // restore stream state
  }
  virtual const char *what(void) const noexcept {
    std::fprintf(stderr, "%s:%ld:%ld: %s\n", fn.c_str(), ln, cn, msg.c_str());
    std::fprintf(stderr, "%5ld | %s\n", ln, line.c_str());
    std::fprintf(stderr, //
                 ("%5s | %" + std::to_string(cn) + "c\n").c_str(), "", '^');
    return "parse error";
  }

  static void check_eof(std::istream &is, const std::string &filename) {
    if (is.eof() || !is.good()) {
      throw parse_error(is, filename, "unexpecteded eof");
    }
  };
};

inline auto read_or_throw(const std::string &filename) {
  if (std::ifstream ifs(filename); ifs) {
    return ifs;
  }
  throw std::runtime_error("cannot read from '" + filename + '\'');
}
}; // namespace __config_detail

struct config : __config_detail::var_t {
  using var_t = __config_detail::var_t;
  using nil_t = __config_detail::nil_t;
  using arr_t = __config_detail::arr_t;
  using obj_t = __config_detail::obj_t;
  using str_t = __config_detail::str_t;
  using dbl_t = __config_detail::dbl_t;
  using int_t = __config_detail::int_t;
  using bol_t = __config_detail::bol_t;
  using tim_t = __config_detail::tim_t;

  using var_t::variant;

  using parse_error = __config_detail::parse_error;

  size_t size(void) const {
    if (auto p = std::get_if<arr_t>(this)) {
      return p->size();
    } else if (auto p = std::get_if<obj_t>(this)) {
      return p->size();
    } else { // throw?
      return 0;
    }
  }

  bool is_arithmetic(void) const {
    using namespace __config_detail;
    return index() == var_i<dbl_t> || index() == var_i<int_t> ||
           index() == var_i<bol_t>;
  }

#define __config_direct_access(t)                                              \
  const t##_t &t(void) const { return std::get<t##_t>(*this); }                \
  t##_t &t(void) {                                                             \
    return std::get<t##_t>(*(is_set() ? this : (*this = t##_t(), this)));      \
  }
  __config_direct_access(arr);
  __config_direct_access(obj);
  __config_direct_access(str);
  double num(std::optional<double> _fallback = {}) const {
    return is_set() ? (double)*this : _fallback.value();
  }

  template <typename T> T &get(void) { return std::get<T>(*this); }
  template <typename T> const T &get(void) const { return std::get<T>(*this); }
  template <size_t I> auto &get(void) { return std::get<I>(*this); }
  template <size_t I> const auto &get(void) const { return std::get<I>(*this); }
#undef __config_direct_access

  config &operator[](size_t i) { return std::get<arr_t>(*this).at(i); }
  const config &operator[](size_t i) const {
    return i < size() ? std::get<arr_t>(*this).at(i) : nil;
  }
  config &operator[](const std::string &path);
  const config &operator[](const std::string &path) const;

  void unset(void) { emplace<__config_detail::nil_t>(); }
  bool is_set(void) const {
    return !std::holds_alternative<__config_detail::nil_t>(*this);
  }
  bool is_deep_set(void) const {
    if (auto pobj = std::get_if<obj_t>(this)) {
      return pobj->end() !=
             std::find_if(pobj->begin(), pobj->end(),
                          [](auto &c) { return c.second.is_deep_set(); });
    } else {
      return !std::holds_alternative<nil_t>(*this);
    }
  }

  template <typename T> bool holds(void) const {
    return std::holds_alternative<T>(*this);
  }

  // convert from initializer_list
  auto &operator=(std::initializer_list<config> l) {
    return l.size() ? *this = arr_t(l.begin(), l.end())
                    : *this = std::monostate();
  }

  // convert to other types
  template <typename T,
            std::enable_if_t<
                !__config_detail::block_implicit<std::remove_cv_t<T>>::value &&
                    !std::is_convertible_v<T, std::string>,
                bool> = true>
  operator T(void) const {
    if constexpr (__config_detail::from_union<T>) {
      if (std::holds_alternative<T>(*this)) {
        return std::get<T>(*this);
      }
    }
    return __config_detail::config_convert_rt<T>(*this);
  }
  // explicitly conversions to avoid ambiguity
  operator std::string(void) const {
    return str(); // TODO: stringify date/time
  }
  explicit operator long(void) const {
    return std::holds_alternative<long>(*this)
               ? std::get<long>(*this)
               : __config_detail::config_convert_rt<long>(*this);
  }

  template <typename T> operator std::vector<T>(void) const {
    return size() ? std::vector<T>(arr().begin(), arr().end())
                  : std::vector<T>();
  }
  template <typename T> operator std::valarray<T>(void) const {
    std::valarray<T> newarr(size());
    if (newarr.size()) {
      auto it = arr().begin();
      for (auto &e : newarr) {
        e = *it++;
      }
    }
    return newarr;
  }

  // set/get with default value
  template <typename T = std::initializer_list<config>>
  auto or_get(const T &default_val) const {
    return is_set() ? *this : config(default_val);
  }
  template <typename T = std::initializer_list<config>>
  auto &or_set(const T &default_val) {
    return is_set() ? *this : *this = default_val;
  }
  template <typename T = std::initializer_list<config>>
  auto &or_set(size_t n, const T &v) {
    return is_set() ? *this : *this = arr_t(n, v);
  }

  // compare operators
  template <typename T> bool operator==(const T &v) const {
    using namespace __config_detail;
    if constexpr (std::is_arithmetic_v<T>) {
      return is_arithmetic() && v == num();
    } else if constexpr (std::is_same_v<T, nil_t>) {
      return index() == var_i<nil_t>;
    } else {
      return index() == var_i<str_t> && v == str();
    }
  }
  template <typename T> friend bool operator==(const T &v, const config &c) {
    return c == v;
  }

  bool operator==(const config &c) const {
    using namespace __config_detail;
    if (this == &c) {
      return true;
    } else if (index() == c.index()) {
      return static_cast<var_t>(*this) == static_cast<var_t>(c);
    } else if (is_arithmetic() && c.is_arithmetic()) {
      return num() == c.num();
    } else {
      return false;
    }
  }

  std::ostream &write_json_inline(std::ostream &os) const;
  std::ostream &write_toml_inline(std::ostream &os) const;
  // std::ostream &write_json(std::ostream &os) const;
  // std::ostream &write_toml(std::ostream &os) const;

  std::istream &parse_json(std::istream &is,
                           const std::string &fn = "<anonymous>");
  std::istream &parse_toml(std::istream &is,
                           const std::string &fn = "<anonymous>");
  std::istream &parse_toml_inline(std::istream &is,
                                  const std::string &fn = "<anonymous>");

  config &parse_json(const std::string &fn) {
    auto ifs = __config_detail::read_or_throw(fn);
    return parse_json(ifs, fn), *this;
  }
  config &parse_toml(const std::string &fn) {
    auto ifs = __config_detail::read_or_throw(fn);
    return parse_toml(ifs, fn), *this;
  }
  config &parse_auto(const std::string &fn) {
    if (fn.substr(fn.size() - 5) == ".json") {
      return parse_json(fn);
    } else if (fn.substr(fn.size() - 5) == ".toml") {
      return parse_toml(fn);
    } else {
      std::vector<parse_error> ve;
      try {
        return parse_json(fn);
      } catch (parse_error &e) {
        e.msg = "when parsed as json: " + e.msg;
        ve.push_back(std::move(e));
      }
      try {
        return parse_json(fn);
      } catch (parse_error &e) {
        e.msg = "when parsed as toml: " + e.msg;
        ve.push_back(std::move(e));
      }
      for (size_t i = 0; i < ve.size() - 1; ++i) {
        ve[i].what();
      }
      throw ve.back();
    }
  }

  const static config nil;
};

inline const config config::nil;

inline config &config::operator[](const std::string &path) {
  using namespace __config_detail;
  auto pcfg = this;
  std::stringstream ss(path);
  while (ss) {
    auto key = read_word(ss, '.');
    ss >> std::ws, ss.peek() == '.' && ss.ignore();
    if (key.size()) {
      if (pcfg->holds<arr_t>()) { // access last element in array
        if (pcfg->arr().empty()) {
          pcfg->arr().emplace_back();
        }
        pcfg = &pcfg->arr().back();
      }
      if (!pcfg->is_set()) {
        *pcfg = obj_t();
      }
      if (pcfg->holds<obj_t>()) {
        for (auto &&[k, v] : pcfg->obj()) {
          if (k == key) {
            pcfg = &v;
            goto __next_component;
          }
        }
        pcfg = &pcfg->obj().emplace_back(key, config()).second;
      }
    }
  __next_component:;
  }
  return *pcfg;
}

inline const config &config::operator[](const std::string &path) const {
  using namespace __config_detail;
  auto pcfg = this;
  std::stringstream ss(path);
  while (ss) {
    auto key = read_word(ss, '.');
    ss >> std::ws, ss.peek() == '.' && ss.ignore();
    if (key.size()) {
      if (pcfg->holds<arr_t>() && pcfg->arr().size()) {
        pcfg = &pcfg->arr().back();
      }
      if (pcfg->holds<obj_t>()) {
        for (auto &&[k, v] : pcfg->obj()) {
          if (k == key) {
            pcfg = &v;
            goto __next_component;
          }
        }
      }
      return nil;
    }
  __next_component:;
  }
  return *pcfg;
}

inline std::ostream &config::write_json_inline(std::ostream &os) const {
  using namespace __config_detail;
  if (std::holds_alternative<nil_t>(*this)) {
    os << "null";
  } else if (auto p = std::get_if<obj_t>(this)) {
    for (size_t i = 0; i < p->size(); ++i) {
      const auto &[k, v] = p->operator[](i);
      os << (i ? ',' : '{');
      write_quoted(os, k);
      os << ':';
      v.write_json_inline(os);
    }
    (p->empty() ? os << '{' : os) << "}";
  } else if (auto p = std::get_if<arr_t>(this)) {
    for (size_t i = 0; i < p->size(); ++i) {
      const auto &e = p->operator[](i);
      os << (i ? ',' : '[');
      e.write_json_inline(os);
    }
    (p->empty() ? os << '[' : os) << "]";
  } else if (auto p = std::get_if<str_t>(this)) {
    write_quoted(os, *p);
  } else if (auto p = std::get_if<dbl_t>(this)) {
    os << *p;
  } else if (auto p = std::get_if<int_t>(this)) {
    os << *p;
  } else if (auto p = std::get_if<bol_t>(this)) {
    os << (*p ? "true" : "false");
  } else {
    throw std::logic_error("impossible branch reached");
  }
  return os;
}

inline std::ostream &config::write_toml_inline(std::ostream &os) const {
  using namespace __config_detail;
  if (std::holds_alternative<nil_t>(*this)) {
    os << "{}"; // XXX toml does not have null type
  } else if (auto p = std::get_if<obj_t>(this)) {
    char sep = '{';
    for (const auto &[k, v] : *p) {
      if (v.is_deep_set()) { // null values are skiped [PERFORMANCE IMPACT?]
        os << sep;
        write_quoted(os, k);
        os << '=';
        v.write_toml_inline(os);
        sep = ',';
      }
    }
    (sep == '{' ? os << '{' : os) << '}';
  } else if (auto p = std::get_if<arr_t>(this)) {
    for (size_t i = 0; i < p->size(); ++i) {
      const auto &e = p->operator[](i);
      os << (i ? ',' : '[');
      e.write_toml_inline(os);
    }
    (p->empty() ? os << '[' : os) << "]";
  } else if (auto p = std::get_if<str_t>(this)) {
    write_quoted(os, *p);
  } else if (auto p = std::get_if<dbl_t>(this)) {
    os << *p;
  } else if (auto p = std::get_if<int_t>(this)) {
    os << *p;
  } else if (auto p = std::get_if<bol_t>(this)) {
    os << (*p ? "true" : "false");
  } else {
    throw std::logic_error("impossible branch reached");
  }
  return os;
}

inline std::istream &config::parse_json(std::istream &is,
                                        const std::string &fn) {
  using namespace __config_detail;
  switch (char c = EOF; is >> c, c) {
  case EOF: {
    return is;
  }
  case '"': {
    auto q = string_unescape(read_quoted(is));
    parse_error::check_eof(is, fn);
    return *this = q, is;
  }
  case '[': {
    arr_t tmp_arr;
    while (is >> std::ws) {
      if (is.peek() == ']' && is.ignore()) {
        break;
      }
      parse_error::check_eof(is, fn);
      tmp_arr.emplace_back().parse_json(is) >> std::ws;
      if (char c = is.get(); c == ']') {
        break;
      } else if (c != ',') {
        is.unget();
        // TODO message
        throw parse_error(is, fn, "expect ',' between json array elements");
      }
    }
    parse_error::check_eof(is, fn);
    return *this = std::move(tmp_arr), is;
  }
  case '{': {
    obj_t tmp_obj;
    while (is >> std::ws) {
      if (is.peek() == '}' && is.ignore()) {
        break;
      }
      if (is.get() != '"') {
        throw parse_error(is, fn, "expect string as object key");
      }
      tmp_obj.emplace_back(string_unescape(read_quoted(is)), config::nil);
      if (is >> std::ws, is.get() != ':') {
        throw parse_error(is, fn, "expect ':' after object key");
      }
      parse_error::check_eof(is, fn);
      tmp_obj.back().second.parse_json(is) >> std::ws;
      if (char c = is.get(); c == '}') {
        break;
      } else if (c != ',') {
        throw parse_error(is, fn, "expect ',' between key-value pairs");
      }
    }
    parse_error::check_eof(is, fn);
    return *this = std::move(tmp_obj), is;
  }
  case 'n': {
    if (is.peek() == 'a' && is.ignore() && is.get() == 'n') {
      return *this = NAN, is;
    } else if (is.get() == 'u' && is.get() == 'l' && is.get() == 'l') {
      return this->unset(), is;
    }
    throw parse_error(is, fn, "unexpected symbol. typo for 'null'?");
  }
  case 't': {
    if (is.get() == 'r' && is.get() == 'u' && is.get() == 'e') {
      return *this = true, is;
    }
    throw parse_error(is, fn, "unexpected symbol. typo for 'true'?");
  }
  case 'f': {
    if (is.get() == 'a' && is.get() == 'l' && //
        is.get() == 's' && is.get() == 'e') {
      return *this = false, is;
    }
    throw parse_error(is, fn, "unexpected symbol. typo for 'false'?");
  }
  default: { // try number
    is.unget();
    dbl_t a;
    is >> a;
    if (a == (int_t)a) {
      *this = (int_t)a;
    } else {
      *this = a;
    }
    return is.clear(), is;
  }
  }
  return is;
}

inline std::istream &config::parse_toml_inline(std::istream &is,
                                               const std::string &fn) {
  using namespace __config_detail;
  const auto _trim = [&is](void) -> auto & { return stream_trim(is, "#"); };
  switch (char c = EOF; _trim() >> c, c) {
  case EOF: {
    return is;
  }
  case '"': { // basic string
    auto pos = is.tellg();
    std::string str;
    if (is.get() == c && is.get() == c) { // multiline """..."""
      char r;
      while (is.get(r)) {
        if (r == c) {
          if (r = is.get(), r == c && is.peek() == c && is.ignore()) {
            while (is.peek() == c && is.ignore()) {
              str += c;
            }
            is.clear();
            break;
          } else {
            str += c, str += r;
          }
        } else if (r == '\\') {
          is.get(r);
          if (r == '\n' || (r == '\r' && is.peek() == '\n')) {
            is >> std::ws;
          } else {
            str += '\\', str += r;
          }
        } else {
          str += r;
        }
      }
      parse_error::check_eof(is, fn);
      // trim newline immediately following the opening delimiter
      r = str[0] == '\r';
      *this = string_unescape(str.data() + (str[r] == '\n' ? r + 1 : 0));
      is.clear();
      return is;
    } else { // non-multiline
      // XXX newline is allowed in our basic string, not by standard
      is.clear(), is.seekg(pos);
      str = read_quoted(is);
      parse_error::check_eof(is, fn);
      return *this = string_unescape(str), is;
    }
  }
  case '\'': { // literal string
    auto pos = is.tellg();
    std::string str;
    if (is.get() == c && is.get() == c) { // multiline '''...'''
      char r;
      while (is.get(r)) {
        if (r == c) {
          if (r = is.get(), r == c && is.peek() == c && is.ignore()) {
            while (is.peek() == c && is.ignore()) {
              str += c;
            }
            is.clear();
            break;
          } else {
            str += c, str += r;
          }
        } else {
          str += r;
        }
      }
      parse_error::check_eof(is, fn);
      // trim newline immediately following the opening delimiter
      r = str[0] == '\r';
      *this = (str.data() + (str[r] == '\n' ? r + 1 : 0));
      return is;
    } else {
      // XXX newline is allowed in our basic string, not by standard
      is.clear(), is.seekg(pos);
      str = read_quoted(is, c, c);
      parse_error::check_eof(is, fn);
      return *this = std::move(str), is;
    }
  }
  case '[': { // array
    arr_t tmp_arr;
    while (_trim()) {
      if (is.peek() == ']' && is.ignore()) {
        break;
      }
      parse_error::check_eof(is, fn);
      tmp_arr.emplace_back(config()).parse_toml_inline(is);
      _trim();
      if (char c = is.get(); c == ']') {
        break;
      } else if (c != ',') {
        is.unget();
        throw parse_error(is, fn, "expect ',' between array elements");
      }
    }
    parse_error::check_eof(is, fn);
    return *this = std::move(tmp_arr), is;
  }
  case '{': { // table
    config tmp_obj = obj_t();
    while (_trim()) {
      if (is.peek() == '}' && is.ignore()) {
        break;
      }
      auto &target = get_context(tmp_obj, is, '=', fn);
      if (is.peek() != '=' || !is.ignore()) {
        throw parse_error(is, fn, "expect '=' after object key");
      }
      _trim(), parse_error::check_eof(is, fn);
      target.parse_toml_inline(is);

      if (_trim().peek() == ',' && is.ignore()) {
        continue;
      } else if (is.peek() == '}' && is.ignore()) {
        break;
      }
      throw parse_error(is, fn, "expect ',' between key-value pairs");
    }
    return *this = std::move(tmp_obj), is;
  }
  case 't': {
    if (is.get() == 'r' && is.get() == 'u' && is.get() == 'e') {
      return *this = true, is;
    }
    throw parse_error(is, fn, "unexpected symbol.");
  }
  case 'f': {
    if (is.get() == 'a' && is.get() == 'l' && //
        is.get() == 's' && is.get() == 'e') {
      return *this = false, is;
    }
    throw parse_error(is, fn, "unexpected symbol.");
  }
  case 'n': {
    if (is.get() == 'a' && is.get() == 'n') {
      return *this = NAN, is;
    }
    throw parse_error(is, fn, "unexpected symbol.");
  }
  case 'i': {
    if (is.get() == 'n' && is.get() == 'f') {
      return *this = HUGE_VAL, is;
    }
    throw parse_error(is, fn, "unexpected symbol.");
  }
  default: {
    int sign = +1;
    if (c == '0') {
      switch (unsigned int nbit = 0; is.peek()) {
        // unsigned integer in different bases
      case 'x': // hex 4bits
        nbit += 1;
        [[fallthrough]];
      case 'o': // oct 3bits
        nbit += 2;
        [[fallthrough]];
      case 'b': // bin 1bits
        nbit += 1;
        is.ignore();
        {
          unsigned long val = 0;
          while (is.get(c)) {
            if (c == '_') {
              continue;
            }
            unsigned int dig = c >= 'a'   ? c - 'a' + 10
                               : c >= 'A' ? c - 'A' + 10
                               : c > '9'  ? -1
                                          : c - '0';
            if (dig >> nbit) { // not a valid char
              is.unget();
              break;
            }
            val = (val << nbit) + dig;
          }
          if (val > std::numeric_limits<int_t>::max()) {
            *this = (dbl_t)val; // warn?
          } else {
            *this = (int_t)val;
          }
          return is.clear(), is;
        }
      default:
        break;
      }
    } else if (::isdigit(c)) {
      is.unget();
    } else if (c == '-' || c == '+') {
      sign = c == '-' ? -1 : +1;
      if (is.peek() == 'n' && is.ignore()) { // +-nan
        if (is.get() == 'a' && is.get() == 'n') {
          return *this = sign * NAN, is;
        }
        throw parse_error(is, fn, "unexpected character.");
      } else if (is.peek() == 'i' && is.ignore()) { // +-inf
        if (is.get() == 'n' && is.get() == 'f') {
          return *this = sign * HUGE_VAL, is;
        }
        throw parse_error(is, fn, "unexpected character.");
      }
    } else {
      throw parse_error(is, fn, "unexpected character.");
    }
    // auto pos = is.tellg();
    uint64_t vint = 0;
    bool in_int = true;
    int vexp = 0;
    int dig_cap = std::numeric_limits<dbl_t>::max_digits10;
    constexpr auto lim = (uint64_t)1 + std::numeric_limits<int_t>::max();
    while (is.get(c)) {
      if (::isdigit(c)) {
        if ((in_int && (vint <= lim)) || (!in_int && dig_cap > 0)) {
          vint *= 10;
          vint += c - '0';
          --dig_cap;
          vexp -= !in_int;
        } else {
          in_int &= dig_cap > 0;
          vexp += in_int;
        }
      } else if (c == '_') { // grouped digits
        continue;
      } else if (in_int && c == '.') {
        in_int = false;
      } else if (c == 'e' || c == 'E') {
        int e = 0;
        is >> e;
        vexp += e;
        in_int = false;
        break;
      } else {
        is.unget();
        break;
      }
    }
    is.clear();
    if (c != ':' && c != '-') {
      if (in_int) {
        if (vint < lim || (vint == lim && sign < 0)) {
          return *this = (int_t)(sign * vint), is;
        }
      }
      // call magical stdlib for consistent floating point error
      std::stringstream ss;
      ss << vint << 'e' << vexp;
      *this = sign * std::atof(ss.str().c_str());
      return is;
    }

    // TODO: date time support after c++20
    { // date/time: current pos: 1970][-01-01
      tim_t t;
      while (is.get(c) && !::isspace(c)) { // by pass for now
        if (!::isalnum(c) && c != '-' && c != '+' && c != '.') {
          break;
        }
      }
      is.unget(), is.clear();
      return *this = t, is;
    }
  }
  }
  return is;
}

inline std::istream &config::parse_toml(std::istream &is,
                                        const std::string &fn) {
  using namespace __config_detail;
  const auto _trim = [&is](void) -> auto & { return stream_trim(is, "#"); };

  config *context = this;

  while (_trim()) {
    if (char c = is.peek(); c == '[' && is.ignore() && _trim()) {
      if (is.peek() == '[' && is.ignore() && _trim()) { // array
        context = &get_context(*this, is, ']', fn);
        if (!context->holds<arr_t>()) {
          *context = arr_t();
        }
        context = &context->arr().emplace_back(obj_t());
        if (_trim().get() != ']' || _trim().get() != ']') {
          throw parse_error(is, fn, "unclosed array bracket");
        }
      } else { // table
        context = &get_context(*this, is, ']', fn);
        if (_trim().get() != ']') {
          throw parse_error(is, fn, "unclosed table bracket");
        }
      }
      parse_error::check_eof(is, fn);
      continue;
    } else if (c == ',') {
      throw parse_error(is, fn, "expect key but get ','");
    } else {
      auto &target = get_context(*context, is, '=', fn);
      if (is.get() != '=') {
        throw parse_error(is, fn, "expect '=' between key-value pair");
      }
      _trim(), parse_error::check_eof(is, fn);
      target.parse_toml_inline(is);
    }
  }

  return is;
}

inline config &__config_detail::get_context(config &root, std::istream &is,
                                            char del, const std::string &fn) {
  std::string path;
  std::getline(stream_trim(is, "#"), path, del);
  if (path.empty()) {
    throw parse_error(is, fn, "expect key");
  }
  is.unget();
  return root[path];
}
}; // namespace yuc
// vi:ft=cpp
