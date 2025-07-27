#include "plotter.h"
#include "supportLib.hpp"

#include <vector>

std::vector<wchar_t>* new_vec_char (const std::wstring& str)
{
    return new std::vector<wchar_t>(str.begin(), str.end());
}

bool GeneratePlotFromFunc(const std::string& filename,
    const std::function<double(double)>& gen, const uint32_t& num, double xmin, double xmax)
{
    std::vector<double> xs;
    std::vector<double> ys;

    double xi = xmin;
    double xstep = (xmax-xmin)/num;

    for (auto i=0; i<num+1; i++)
    {
        xs.push_back(xi);
        ys.push_back(gen(xs.back()));
        xi+=xstep;
    }

    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = new std::vector<double>(xs);
    series->ys = new std::vector<double>(ys);
    series->linearInterpolation = true;
    series->lineType = new_vec_char(plot_data.line_type);
    series->lineThickness = 2;
    series->color = CreateRGBColor(plot_data.rgb[0], plot_data.rgb[1], plot_data.rgb[2]);

    return GeneratePlot(filename, series);
}

bool GeneratePlotFromPoints(const std::string& filename, 
                            const std::vector<double>& xs, 
                            const std::vector<double>& ys)
{
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = new std::vector<double>(xs);
    series->ys = new std::vector<double>(ys);
    series->linearInterpolation = true;
    series->lineType = new_vec_char(plot_data.line_type);
    series->lineThickness = 2;
    series->color = CreateRGBColor(plot_data.rgb[0], plot_data.rgb[1], plot_data.rgb[2]);

    return GeneratePlot(filename, series);
}

bool CalculateBounds(ScatterPlotSettings* settings, ScatterPlotSeries* series)
{
    if(!settings || !series) 
        return false;

    settings->xMin = FLT_MAX;
    settings->xMax = -FLT_MAX;
    settings->yMin = FLT_MAX;
    settings->yMax = -FLT_MAX;

    for (int i = 0; i < series->xs->size(); i++)
    {
        settings->xMin = std::min((*series->xs)[i], settings->xMin);
        settings->xMax = std::max((*series->xs)[i], settings->xMax);
        settings->yMin = std::min((*series->ys)[i], settings->yMin);
        settings->yMax = std::max((*series->ys)[i], settings->yMax);
    }

    plot_data.range_x_min = settings->xMin;
    plot_data.range_x_max = settings->xMax;
    plot_data.range_y_min = settings->yMin;
    plot_data.range_y_max = settings->yMax;
    return true;
}

bool GeneratePlot(const std::string& filename, ScatterPlotSeries *series) {

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();

    CalculateBounds(settings, series);

    settings->width = plot_data.pix_x;
    settings->height = plot_data.pix_y;
    settings->autoBoundaries = false;
    settings->autoPadding = false;
    settings->xPadding=plot_data.pad_x;
    settings->yPadding=plot_data.pad_y;
    settings->title = new_vec_char(plot_data.plot_name);
    settings->xLabel = new_vec_char(L"X axis");
    settings->yLabel = new_vec_char(L"Y axis");
    settings->scatterPlotSeries = new std::vector<ScatterPlotSeries*> {series};

    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = new StringReference();
    bool success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

    if (success)
    {
        std::vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, filename);
        DeleteImage(imageReference->image);
    }

    return success;
}

bool GenerateEmptyPlot(const std::string& filename) {

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();

    settings->xMin = plot_data.range_x_min;
    settings->xMax = plot_data.range_x_max;
    settings->yMin = plot_data.range_y_min;
    settings->yMax = plot_data.range_y_max;

    settings->width = plot_data.pix_x;
    settings->height = plot_data.pix_y;
    settings->autoBoundaries = false;
    settings->autoPadding = false;
    settings->xPadding=plot_data.pad_x;
    settings->yPadding=plot_data.pad_y;
    settings->title = new_vec_char(L"Empty");
    settings->xLabel = new_vec_char(L"X axis");
    settings->yLabel = new_vec_char(L"Y axis");

    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = new StringReference();
    bool success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

    if (success)
    {
        std::vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, filename);
        DeleteImage(imageReference->image);
    }

    return success;
}

bool GenerateSimplePlot(const std::string& filename, std::vector<double>& xs,
                        std::vector<double>& ys)
{

    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

    StringReference* errorMessage = new StringReference();

    bool success = DrawScatterPlot(imageReference, plot_data.pix_x,
                                   plot_data.pix_y, &xs, &ys, errorMessage);

    if (success)
    {
        std::vector<double>* pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, filename);
        DeleteImage(imageReference->image);
    }

    return success;
}