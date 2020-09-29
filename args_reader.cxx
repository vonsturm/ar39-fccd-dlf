// argument reader
//
// Author  : Katharina von Sturm
// Date    : 29.09.2020
// Usage   : fetch_arg(<arguments>, <identifier>, <var>)
//  searches <arguments> for <identifier> and reads the next element into <var>
// Exception : if <var> is of bool type and <identifier> is found <var> is set to true
// Return value : returns `true` if <identifier> was found, `false` otherwise

#include <vector>
#include <algorithm>

template<typename T>
bool fetch_arg(const std::vector<std::string> & args, std::string opt, T & var) {
  auto result = find(args.begin(), args.end(), opt);
  if (result != args.end()) {
         if constexpr (std::is_same<T, int>        ::value) var = stoi(*(result+1));
    else if constexpr (std::is_same<T, float>      ::value) var = stof(*(result+1));
    else if constexpr (std::is_same<T, double>     ::value) var = stod(*(result+1));
    else if constexpr (std::is_same<T, std::string>::value) var = *(result+1);
    else if constexpr (std::is_same<T, bool>       ::value) var = true;
    return true;
  }
  return false;
}