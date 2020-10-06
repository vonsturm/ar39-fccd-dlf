// The MIT License (MIT)
//
// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// ============================================================================
//    ___   ___   __    __      ___   ____   ___   ___   ____  ___
//   / /_\ | |_) / /`_ ( (` __ | |_) | |_   / /_\ | | \ | |_  | |_)
//  /_/  \ |_| \ \_\_/ _)_)    |_| \ |_|__ /_/  \ |_|_/ |_|__ |_| \
//
// Author  : Katharina von Sturm
// Date    : 29.09.2020
//
// Usage   : fetch_arg(<arguments>, <identifier>, <var>)
//    searches <arguments> for <identifier> and reads the next element into <var>
//    if <var> is a vector the number of elements read is equal to the vector size
//    in this case only same_type arguments can be stored into the vector
//
// Exception : if <var> is of bool type and <identifier> is found <var> is set to true
//
// Return value : returns `true` if <identifier> was found, `false` otherwise
//
// ============================================================================

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <type_traits>

template<typename T>
bool fetch_arg(const std::vector<std::string> & args, std::string identifier, T & var) noexcept {
  auto result = find(args.begin(), args.end(), identifier);
  if (result != args.end()) {
         if constexpr (std::is_same<T, int>         ::value) var = stoi(*(result+1));
    else if constexpr (std::is_same<T, float>       ::value) var = stof(*(result+1));
    else if constexpr (std::is_same<T, double>      ::value) var = stod(*(result+1));
    else if constexpr (std::is_same<T, std::string> ::value) var = *(result+1);
    else if constexpr (std::is_same<T, bool>        ::value) var = true;
    return true;
  }
  return false;
}

template<typename T>
bool fetch_arg(const std::vector<std::string> & args, std::string identifier, std::vector<T> & var) noexcept {
  auto result = find(args.begin(), args.end(), identifier);
  if (result != args.end()) {
      int i = 0;
      for (auto & v : var) {
             if constexpr (std::is_same<T, int>         ::value) v = stoi(*(result+1+i++));
        else if constexpr (std::is_same<T, float>       ::value) v = stof(*(result+1+i++));
        else if constexpr (std::is_same<T, double>      ::value) v = stod(*(result+1+i++));
        else if constexpr (std::is_same<T, std::string> ::value) v = *(result+1+i++);
        else if constexpr (std::is_same<T, bool>        ::value) v = true;
      }
    return true;
  }
  return false;
}
