#include <sstream>
#include "plot.h"

Plot::Plot(std::string name,
           std::string x_label,
           std::string y_label,
           std::vector<int> xs) :
    name_(std::move(name)),
    x_label_(std::move(x_label)),
    y_label_(std::move(y_label)),
    xs_(std::move(xs)) {}

void Plot::AddPlotLine(PlotLine plot_line) {
  lines_.push_back(std::move(plot_line));
}

std::string Plot::ToString() {
  std::stringstream ss;
  ss << name_ << '\n' << x_label_ << '\n' << y_label_ << '\n';
  for (auto x: xs_) {
    ss << x << ' ';
  }
  ss << '\n';
  for (const auto& line: lines_) {
    ss << line.GetName() << '\n';
    for (auto x: xs_) {
      if (line.GetValues().count(x) != 0) {
        ss << line.GetValues().at(x);
      } else {
        ss << -1;
      }
      ss << ' ';
    }
    ss << '\n';
  }
  ss << "\n===\n";
  return ss.str();
}
