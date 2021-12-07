#pragma once

#include <string>
#include <vector>
#include <map>

class PlotLine {
 public:
  explicit PlotLine(std::string name);
  std::map<int, double>& GetValues();
  const std::map<int, double>& GetValues() const;
  const std::string& GetName() const;
  void AddValue(int x, double y);

 private:
  std::string name_;
  std::map<int, double> values_;
};
