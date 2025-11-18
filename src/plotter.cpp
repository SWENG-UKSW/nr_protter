#include "plotter.h"
#include "supportLib.hpp"

#include <cfloat>
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

bool GenerateContinuousPlotFromFunc(const std::string& filename,
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

    return ContinuousPlot(filename, series);
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


namespace {

struct Rectangle{
    double x1;
    double x2;
    double y1;
    double y2;
};

void CompBoundariesBasedOnSettings(ScatterPlotSettings *settings, Rectangle *boundaries){
    ScatterPlotSeries *sp;
    double plot, xMin, xMax, yMin, yMax;

    if(settings->scatterPlotSeries->size() >= 1.0){
        xMin = GetMinimum(settings->scatterPlotSeries->at(0)->xs);
        xMax = GetMaximum(settings->scatterPlotSeries->at(0)->xs);
        yMin = GetMinimum(settings->scatterPlotSeries->at(0)->ys);
        yMax = GetMaximum(settings->scatterPlotSeries->at(0)->ys);
    }else{
        xMin =  -10.0;
        xMax = 10.0;
        yMin =  -10.0;
        yMax = 10.0;
    }

    if( !settings->autoBoundaries ){
        xMin = settings->xMin;
        xMax = settings->xMax;
        yMin = settings->yMin;
        yMax = settings->yMax;
    }else{
        for(plot = 1.0; plot < settings->scatterPlotSeries->size(); plot = plot + 1.0){
            sp = settings->scatterPlotSeries->at(plot);

            xMin = fmin(xMin, GetMinimum(sp->xs));
            xMax = fmax(xMax, GetMaximum(sp->xs));
            yMin = fmin(yMin, GetMinimum(sp->ys));
            yMax = fmax(yMax, GetMaximum(sp->ys));
        }
    }

    boundaries->x1 = xMin;
    boundaries->y1 = yMin;
    boundaries->x2 = xMax;
    boundaries->y2 = yMax;
}

bool AmendScatterPlotFromSettings(RGBABitmapImageReference *canvasReference, ScatterPlotSettings *settings, StringReference *errorMessage){
    double xMin, xMax, yMin, yMax, xLength, yLength, i, x, y, xPrev, yPrev, px, py, pxPrev, pyPrev, originX, originY, p, l, plot;
    Rectangle boundaries;
    double xPadding, yPadding, originXPixels, originYPixels;
    double xPixelMin, yPixelMin, xPixelMax, yPixelMax, xLengthPixels, yLengthPixels, axisLabelPadding;
    NumberReference nextRectangle, x1Ref, y1Ref, x2Ref, y2Ref, patternOffset;
    bool prevSet, success;
    RGBA *gridLabelColor;
    RGBABitmapImage *canvas;
    std::vector<double> *xs, *ys;
    bool linearInterpolation;
    ScatterPlotSeries *sp;
    std::vector<double> *xGridPositions, *yGridPositions;
    StringArrayReference *xLabels, *yLabels;
    NumberArrayReference *xLabelPriorities, *yLabelPriorities;
    std::vector<Rectangle*> *occupied;
    std::vector<bool> *linePattern;
    bool originXInside, originYInside, textOnLeft, textOnBottom;
    double originTextX, originTextY, originTextXPixels, originTextYPixels, side;

    canvas = canvasReference->image;
    patternOffset.numberValue=0.0;

    success = ScatterPlotFromSettingsValid(settings, errorMessage);

    if(success){


        CompBoundariesBasedOnSettings(settings, &boundaries);
        xMin = boundaries.x1;
        yMin = boundaries.y1;
        xMax = boundaries.x2;
        yMax = boundaries.y2;

        /* If zero, set to defaults. */
        if(xMin - xMax == 0.0){
            xMin = 0.0;
            xMax = 10.0;
        }

        if(yMin - yMax == 0.0){
            yMin = 0.0;
            yMax = 10.0;
        }

        xLength = xMax - xMin;
        yLength = yMax - yMin;

        if(settings->autoPadding){
            xPadding = floor(GetDefaultPaddingPercentage()*settings->width);
            yPadding = floor(GetDefaultPaddingPercentage()*settings->height);
        }else{
            xPadding = settings->xPadding;
            yPadding = settings->yPadding;
        }

        /* Compute origin information. */
        originYInside = yMin < 0.0 && yMax > 0.0;
        originY = 0.0;
        if(settings->xAxisAuto){
            if(originYInside){
                originY = 0.0;
            }else{
                originY = yMin;
            }
        }else{
            if(settings->xAxisTop){
                originY = yMax;
            }
            if(settings->xAxisBottom){
                originY = yMin;
            }
        }
        xPixelMin = xPadding;
        yPixelMin = yPadding;
        xPixelMax = settings->width - xPadding;
        yPixelMax = settings->height - yPadding;
        xLengthPixels = xPixelMax - xPixelMin;
        yLengthPixels = yPixelMax - yPixelMin;

        originYPixels = MapYCoordinate(originY, yMin, yMax, yPixelMin, yPixelMax);

        originXInside = xMin < 0.0 && xMax > 0.0;
        originX = 0.0;
        if(settings->yAxisAuto){
            if(originXInside){
                originX = 0.0;
            }else{
                originX = xMin;
            }
        }else{
            if(settings->yAxisLeft){
                originX = xMin;
            }
            if(settings->yAxisRight){
                originX = xMax;
            }
        }
        originXPixels = MapXCoordinate(originX, xMin, xMax, xPixelMin, xPixelMax);

        if(originYInside){
            originTextY = 0.0;
        }else{
            originTextY = yMin + yLength/2.0;
        }
        originTextYPixels = MapYCoordinate(originTextY, yMin, yMax, yPixelMin, yPixelMax);

        if(originXInside){
            originTextX = 0.0;
        }else{
            originTextX = xMin + xLength/2.0;
        }
        originTextXPixels = MapXCoordinate(originTextX, xMin, xMax, xPixelMin, xPixelMax);

        /* Draw points */
        for(plot = 0.0; plot < settings->scatterPlotSeries->size(); plot = plot + 1.0){
            sp = settings->scatterPlotSeries->at(plot);

            xs = sp->xs;
            ys = sp->ys;
            linearInterpolation = sp->linearInterpolation;

            if(linearInterpolation){
                prevSet = false;
                xPrev = 0.0;
                yPrev = 0.0;
                for(i = 0.0; i < xs->size(); i = i + 1.0){
                    x = xs->at(i);
                    y = ys->at(i);

                    if(prevSet){
                        x1Ref.numberValue = xPrev;
                        y1Ref.numberValue = yPrev;
                        x2Ref.numberValue = x;
                        y2Ref.numberValue = y;

                        success = CropLineWithinBoundary(&x1Ref, &y1Ref, 
                            &x2Ref, &y2Ref, xMin, xMax, yMin, yMax);

                        if(success){
                            pxPrev = floor(MapXCoordinate(x1Ref.numberValue, xMin, xMax, xPixelMin, xPixelMax));
                            pyPrev = floor(MapYCoordinate(y1Ref.numberValue, yMin, yMax, yPixelMin, yPixelMax));
                            px = floor(MapXCoordinate(x2Ref.numberValue, xMin, xMax, xPixelMin, xPixelMax));
                            py = floor(MapYCoordinate(y2Ref.numberValue, yMin, yMax, yPixelMin, yPixelMax));

                            if(aStringsEqual(sp->lineType, toVector(L"solid")) && sp->lineThickness == 1.0){
                                DrawLine1px(canvas, pxPrev, pyPrev, px, py, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"solid"))){
                                DrawLine(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"dashed"))){
                                linePattern = GetLinePattern1();
                                DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, &patternOffset, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"dotted"))){
                                linePattern = GetLinePattern2();
                                DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, &patternOffset, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"dotdash"))){
                                linePattern = GetLinePattern3();
                                DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, &patternOffset, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"longdash"))){
                                linePattern = GetLinePattern4();
                                DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, &patternOffset, sp->color);
                            }else if(aStringsEqual(sp->lineType, toVector(L"twodash"))){
                                linePattern = GetLinePattern5();
                                DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, &patternOffset, sp->color);
                            }
                        }
                    }

                    prevSet = true;
                    xPrev = x;
                    yPrev = y;
                }
            }else{
                for(i = 0.0; i < xs->size(); i = i + 1.0){
                    x = xs->at(i);
                    y = ys->at(i);

                    if(x > xMin && x < xMax && y > yMin && y < yMax){

                        x = floor(MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax));
                        y = floor(MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax));

                        if(aStringsEqual(sp->pointType, toVector(L"crosses"))){
                            DrawPixel(canvas, x, y, sp->color);
                            DrawPixel(canvas, x + 1.0, y, sp->color);
                            DrawPixel(canvas, x + 2.0, y, sp->color);
                            DrawPixel(canvas, x - 1.0, y, sp->color);
                            DrawPixel(canvas, x - 2.0, y, sp->color);
                            DrawPixel(canvas, x, y + 1.0, sp->color);
                            DrawPixel(canvas, x, y + 2.0, sp->color);
                            DrawPixel(canvas, x, y - 1.0, sp->color);
                            DrawPixel(canvas, x, y - 2.0, sp->color);
                        }else if(aStringsEqual(sp->pointType, toVector(L"circles"))){
                            DrawCircle(canvas, x, y, 3.0, sp->color);
                        }else if(aStringsEqual(sp->pointType, toVector(L"dots"))){
                            DrawFilledCircle(canvas, x, y, 3.0, sp->color);
                        }else if(aStringsEqual(sp->pointType, toVector(L"triangles"))){
                            DrawTriangle(canvas, x, y, 3.0, sp->color);
                        }else if(aStringsEqual(sp->pointType, toVector(L"filled triangles"))){
                            DrawFilledTriangle(canvas, x, y, 3.0, sp->color);
                        }else if(aStringsEqual(sp->pointType, toVector(L"pixels"))){
                            DrawPixel(canvas, x, y, sp->color);
                        }
                    }
                }
            }
        }
    }

    return success;
}

}

static RGBABitmapImageReference* gimageReference=nullptr;
bool ContinuousPlot(const std::string& filename, ScatterPlotSeries *series) {

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
    settings->title = new_vec_char(plot_data.plot_name);
    settings->xLabel = new_vec_char(L"X axis");
    settings->yLabel = new_vec_char(L"Y axis");
    settings->scatterPlotSeries = new std::vector<ScatterPlotSeries*> {series};

    bool firstInLine=false;
    if(!gimageReference)
    {
        gimageReference = CreateRGBABitmapImageReference();
        gimageReference->image = CreateImage(settings->width, settings->height, GetWhite());

        firstInLine=true;
    }

    StringReference *errorMessage = new StringReference();
    bool success = false;

    if (firstInLine)
    {
        success = DrawScatterPlotFromSettings(gimageReference, settings, errorMessage);
    }
    else {
        success = AmendScatterPlotFromSettings(gimageReference, settings, errorMessage);
    }



    if (success)
    {
        std::vector<double> *pngdata = ConvertToPNG(gimageReference->image);
        WriteToFile(pngdata, filename);
    }

    return success;
}

void FinishContinuousPlot()
{
    DeleteImage(gimageReference->image);
    gimageReference=nullptr;
}
