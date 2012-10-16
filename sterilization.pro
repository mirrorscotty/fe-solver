######################################################################
# Automatically generated by qmake (2.01a) Mon Aug 27 02:21:37 2012
######################################################################

QMAKEFEATURES += /usr/lib64/qt4/features6
CONFIG += qwt
unix:LIBS += -lqwt6
unix:INCLUDEPATH += /usr/include/qwt6
win32:LIBS += C:/Qwt-6.0.1/lib/qwt.dll
win32:INCLUDEPATH += C:/Qwt-6.0.1/include

QMAKE_CXXFLAGS = -ggdb

TEMPLATE = app
TARGET = solver
OBJECTS_DIR = ./tmp

DEPENDPATH += . \
              scaling \
              fem \
              gui/heating \
              matrix \
              mesh \
              problems \
              extras/doc \
              material-data/choi-okos
INCLUDEPATH += . \
               matrix \
               scaling \
               fem \
               mesh \
               gui \
               material-data/choi-okos

# Input
HEADERS += arclength.h \
           basis.h \
           integrate.h \
           isoparam.h \
           fem/solve.h \
           output.h \
           scaling/scaling_ht.h \
           fem/auxsoln.h \
           fem/finite-element.h \
           fem/finite-element1d.h \
           fem/solution.h \
           gui/heating/solver.h \
           gui/heating/heat-gui.h \
           matrix/matrix.h \
           mesh/mesh1d.h \
           mesh/mesh2d.h \
           material-data/choi-okos/choi-okos.h \
           material-data/choi-okos/datafile.h \
           matrix/mtxsolver.h \
           matrix/bandmatrix.h \
           about.h
FORMS += gui/heating/solver.ui
SOURCES += basis.c \
           integrate.c \
           isoparam.c \
           output.c \
           fem/solve.c \
           scaling/scaling_ht.c \
           fem/auxsoln.c \
           fem/finite-element.c \
           fem/finite-element1d.c \
           fem/solution.c \
           gui/heating/main.cpp \
           gui/heating/solver.cpp \
           gui/heating/heat-gui.c \
           matrix/matrix.c \
           matrix/vector.c \
           mesh/mesh1d.c \
           mesh/mesh2d.c \
           material-data/choi-okos/choi-okos.c \
           material-data/choi-okos/datafile.c \
           matrix/mtxsolver.c \
           matrix/bandmatrix.c \
           about.c
