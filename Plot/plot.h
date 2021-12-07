#pragma once

#include "plot_line.h"

class Plot {
 public:
  Plot(std::string name,
       std::string x_label,
       std::string y_label,
       std::vector<int> xs);

  void AddPlotLine(PlotLine plot_line);

  std::string ToString();

 private:
  std::string name_;
  std::string x_label_;
  std::string y_label_;
  std::vector<int> xs_;
  std::vector<PlotLine> lines_;
};
