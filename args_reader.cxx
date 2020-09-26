#include <vector>
#include <algorithm>

template<typename T>
void fetch_arg(const std::vector<std::string> & args, std::string opt, T & var) {
  auto result = find(args.begin(), args.end(), opt);
  if (result != args.end()) {
         if constexpr (std::is_same<T, int>        ::value) var = stoi(*(result+1));
    else if constexpr (std::is_same<T, float>      ::value) var = stof(*(result+1));
    else if constexpr (std::is_same<T, double>     ::value) var = stod(*(result+1));
    else if constexpr (std::is_same<T, std::string>::value) var = *(result+1);
    else if constexpr (std::is_same<T, bool>       ::value) var = true;
  }
}