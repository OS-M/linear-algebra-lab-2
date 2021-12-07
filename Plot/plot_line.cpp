#include "plot_line.h"

PlotLine::PlotLine(std::string name) : name_(std::move(name)) {}

std::map<int, double>& PlotLine::GetValues() {
  return values_;
}

const std::map<int, double>& PlotLine::GetValues() const {
  return values_;
}
const std::string& PlotLine::GetName() const {
  return name_;
}

void PlotLine::AddValue(int x, double y) {
  values_[x] = y;
}
