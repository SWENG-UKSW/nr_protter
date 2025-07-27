#ifndef PLOTTER_H
#define PLOTTER_H

#include <functional>
#include <string>
#include "pbPlots.hpp"

struct PlotData
{
    int pix_x=0;
    int pix_y=0;

    int pad_x=0;
    int pad_y=0;

    float range_x_min=-10.f;
    float range_x_max=10.f;

    float range_y_min=-10.f;
    float range_y_max=10.f;

    std::vector<double> xs;
    std::vector<double> ys;

    std::wstring line_type=L"dashed";
    std::wstring plot_name=L"nothing";

    double rgb[3]={1.0, 1.0, 1.0};
};

extern PlotData plot_data;

bool GeneratePlotFromFunc(const std::string& filename, 
    const std::function<double(double)> & gen, const uint32_t & num, double xmin, double xmax);
    bool GeneratePlotFromPoints(const std::string& filename, 
    const std::vector<double> &xs, const std::vector<double> &ys);
bool GeneratePlot(const std::string& filename, ScatterPlotSeries *series);
bool GenerateSimplePlot(const std::string& filename, 
    std::vector<double>& xs, 
    std::vector<double>& ys);
bool GenerateEmptyPlot(const std::string& filename);


#endif // PLOTTER_H