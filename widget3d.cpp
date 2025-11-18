#include "widget3d.h"

network *widget3d::getNet() const
{
    return net;
}

void widget3d::setNet(network *value)
{
    net = value;
}

bool widget3d::getLight() const
{
    return light;
}

void widget3d::setLight(bool value)
{
    light = value;
}
bool widget3d::getPoreBodies() const
{
    return poreBodies;
}

void widget3d::setPoreBodies(bool value)
{
    poreBodies = value;
}
bool widget3d::getNodeBodies() const
{
    return nodeBodies;
}

void widget3d::setNodeBodies(bool value)
{
    nodeBodies = value;
}
bool widget3d::getOilVisible() const
{
    return oilVisible;
}

void widget3d::setOilVisible(bool value)
{
    oilVisible = value;
}
bool widget3d::getWaterVisible() const
{
    return waterVisible;
}

void widget3d::setWaterVisible(bool value)
{
    waterVisible = value;
}
bool widget3d::getGasVisible() const
{
    return gasVisible;
}

void widget3d::setGasVisible(bool value)
{
    gasVisible = value;
}
bool widget3d::getWaterWetVisible() const
{
    return waterWetVisible;
}

void widget3d::setWaterWetVisible(bool value)
{
    waterWetVisible = value;
}
bool widget3d::getOilWetVisible() const
{
    return oilWetVisible;
}

void widget3d::setOilWetVisible(bool value)
{
    oilWetVisible = value;
}
double widget3d::getXRot() const
{
    return xRot;
}

void widget3d::setXRot(double value)
{
    xRot = value;
}
double widget3d::getYRot() const
{
    return yRot;
}

void widget3d::setYRot(double value)
{
    yRot = value;
}
double widget3d::getZRot() const
{
    return zRot;
}

void widget3d::setZRot(double value)
{
    zRot = value;
}
double widget3d::getXTran() const
{
    return xTran;
}

void widget3d::setXTran(double value)
{
    xTran = value;
}
double widget3d::getYTran() const
{
    return yTran;
}

void widget3d::setYTran(double value)
{
    yTran = value;
}
double widget3d::getScale() const
{
    return scale;
}

void widget3d::setScale(double value)
{
    scale = value;
}
double widget3d::getAspect() const
{
    return aspect;
}

void widget3d::setAspect(double value)
{
    aspect = value;
}

bool widget3d::getAnimation() const
{
    return animation;
}

void widget3d::setAnimation(bool value)
{
    animation = value;
}
double widget3d::getXInitTran() const
{
    return xInitTran;
}

void widget3d::setXInitTran(double value)
{
    xInitTran = value;
}
double widget3d::getYInitTran() const
{
    return yInitTran;
}

void widget3d::setYInitTran(double value)
{
    yInitTran = value;
}
double widget3d::getZInitTran() const
{
    return zInitTran;
}

void widget3d::setZInitTran(double value)
{
    zInitTran = value;
}
double widget3d::getYInitRot() const
{
    return yInitRot;
}

void widget3d::setYInitRot(double value)
{
    yInitRot = value;
}
double widget3d::getXInitRot() const
{
    return xInitRot;
}

void widget3d::setXInitRot(double value)
{
    xInitRot = value;
}
bool widget3d::getCutX() const
{
    return cutX;
}

void widget3d::setCutX(bool value)
{
    cutX = value;
}
bool widget3d::getCutY() const
{
    return cutY;
}

void widget3d::setCutY(bool value)
{
    cutY = value;
}
bool widget3d::getCutZ() const
{
    return cutZ;
}

void widget3d::setCutZ(bool value)
{
    cutZ = value;
}
double widget3d::getCutXValue() const
{
    return cutXValue;
}

void widget3d::setCutXValue(double value)
{
    cutXValue = value;
}
double widget3d::getCutYValue() const
{
    return cutYValue;
}

void widget3d::setCutYValue(double value)
{
    cutYValue = value;
}
double widget3d::getCutZValue() const
{
    return cutZValue;
}

void widget3d::setCutZValue(double value)
{
    cutZValue = value;
}

std::string widget3d::getPhasePoresPath() const
{
    return phasePoresPath;
}

void widget3d::setPhasePoresPath(const std::string &value)
{
    phasePoresPath = value;
}
std::string widget3d::getPhaseNodesPath() const
{
    return phaseNodesPath;
}

void widget3d::setPhaseNodesPath(const std::string &value)
{
    phaseNodesPath = value;
}
bool widget3d::getRender() const
{
    return render;
}

void widget3d::setRender(bool value)
{
    render = value;
}
bool widget3d::getLoad() const
{
    return load;
}

void widget3d::setLoad(bool value)
{
    load = value;
}
int widget3d::getTotalImages() const
{
    return totalImages;
}

void widget3d::setTotalImages(int value)
{
    totalImages = value;
}
int widget3d::getImageCount() const
{
    return imageCount;
}

void widget3d::setImageCount(int value)
{
    imageCount = value;
}







































bool widget3d::getTissue() const
{
    return tissue;
}

void widget3d::setTissue(bool value)
{
    tissue = value;
}
bool widget3d::getCutXT() const
{
    return cutXT;
}

void widget3d::setCutXT(bool value)
{
    cutXT = value;
}
bool widget3d::getCutYT() const
{
    return cutYT;
}

void widget3d::setCutYT(bool value)
{
    cutYT = value;
}
bool widget3d::getCutZT() const
{
    return cutZT;
}

void widget3d::setCutZT(bool value)
{
    cutZT = value;
}
double widget3d::getCutXTValue() const
{
    return cutXTValue;
}

void widget3d::setCutXTValue(double value)
{
    cutXTValue = value;
}
double widget3d::getCutYTValue() const
{
    return cutYTValue;
}

void widget3d::setCutYTValue(double value)
{
    cutYTValue = value;
}
double widget3d::getCutZTValue() const
{
    return cutZTValue;
}

void widget3d::setCutZTValue(double value)
{
    cutZTValue = value;
}
std::string widget3d::getPhaseBlockPath() const
{
    return phaseBlockPath;
}

void widget3d::setPhaseBlockPath(const std::string &value)
{
    phaseBlockPath = value;
}
bool widget3d::getParticles() const
{
    return particles;
}

void widget3d::setParticles(bool value)
{
    particles = value;
}

double widget3d::getCutXValue2() const
{
    return cutXValue2;
}

void widget3d::setCutXValue2(double value)
{
    cutXValue2 = value;
}

double widget3d::getCutYValue2() const
{
    return cutYValue2;
}

void widget3d::setCutYValue2(double value)
{
    cutYValue2 = value;
}

double widget3d::getCutZValue2() const
{
    return cutZValue2;
}

void widget3d::setCutZValue2(double value)
{
    cutZValue2 = value;
}

double widget3d::getCutXTValue2() const
{
    return cutXTValue2;
}

void widget3d::setCutXTValue2(double value)
{
    cutXTValue2 = value;
}

double widget3d::getCutYTValue2() const
{
    return cutYTValue2;
}

void widget3d::setCutYTValue2(double value)
{
    cutYTValue2 = value;
}

double widget3d::getCutZTValue2() const
{
    return cutZTValue2;
}

void widget3d::setCutZTValue2(double value)
{
    cutZTValue2 = value;
}

bool widget3d::getDrugView() const
{
    return drugView;
}

void widget3d::setDrugView(bool value)
{
    drugView = value;
}

bool widget3d::getOxygenView() const
{
    return oxygenView;
}

void widget3d::setOxygenView(bool value)
{
    oxygenView = value;
}
bool widget3d::getAxis() const
{
    return axis;
}

void widget3d::setAxis(bool value)
{
    axis = value;
}
bool widget3d::getCells() const
{
    return cells;
}

void widget3d::setCells(bool value)
{
    cells = value;
}
bool widget3d::getTAFView() const
{
    return TAFView;
}

void widget3d::setTAFView(bool value)
{
    TAFView = value;
}
bool widget3d::getFNView() const
{
    return FNView;
}

void widget3d::setFNView(bool value)
{
    FNView = value;
}
bool widget3d::getMDEView() const
{
    return MDEView;
}

void widget3d::setMDEView(bool value)
{
    MDEView = value;
}













