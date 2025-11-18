#-------------------------------------------------
#
# Project created by QtCreator 2014-01-02T17:45:02
#
#-------------------------------------------------

QT       += core

QT       += gui

QT       += opengl

QT       += printsupport

TARGET = numPTI
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++11

TEMPLATE = app

win32 {
    LIBS += -lopengl32 $$PWD/Glew/glew32.dll
    LIBS += -L$$PWD/Glew/ -lglew32
    LIBS += -L$$PWD/Glew/ -lglew32s
}

unix {
    LIBS += -lGLEW
    QMAKE_CXXFLAGS += -std=c++11
}

SOURCES += main.cpp \
    pore.cpp \
    node.cpp \
    network.cpp \
    mainwindow.cpp \
    cluster.cpp \
    worker.cpp \
    widget3d.cpp \
    experimental.cpp \
    hoshenKopelmann.cpp \
    loadData.cpp \
    element.cpp \
    solver.cpp \
    generationRegular.cpp \
    misc.cpp \
    qcustomplot.cpp \
    block.cpp \
    particle.cpp \
    artificial.cpp \
    drugflow.cpp \
    particleflow.cpp \
    angiogenesis.cpp \
    Cell.cpp \
    Param.cpp \
    PDE.cpp \
    Mesh.cpp \
    Fibre.cpp \
    cellProliferation.cpp \
    angioFlow.cpp \
    parentvessel.cpp \
    retina.cpp

HEADERS += \
    pore.h \
    node.h \
    network.h \
    tools.h \
    mainwindow.h \
    widget3d.h \
    cluster.h \
    worker.h \
    element.h \
    qcustomplot.h \
    block.h \
    particle.h \
    Cell.h \
    Param.h \
    PDE.h \
    Mesh.h \
    Fibre.h

OTHER_FILES +=

FORMS += \
    mainwindow.ui

CONFIG += warn_off

RESOURCES += \
    resource.qrc
