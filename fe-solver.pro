######################################################################
# Automatically generated by qmake (2.01a) Mon Aug 27 02:21:37 2012
######################################################################

QMAKEFEATURES += /usr/lib64/qt4/features6
CONFIG += qwt debug
unix:LIBS += -lqwt6
unix:INCLUDEPATH += /usr/include/qwt6
win32:LIBS += C:/Qwt-6.0.1/lib/qwt.dll
win32:INCLUDEPATH += C:/Qwt-6.0.1/include

QMAKE_CXXFLAGS = -ggdb

TEMPLATE = app
TARGET = solver
OBJECTS_DIR = ./tmp

DEPENDPATH += . \
#              extras \
              scaling \
              fem \
              gui \
              matrix \
              mesh \
              problems \
              extras/doc \
              material-data/can
#              material-data/drying \
#              material-data/freezing
INCLUDEPATH += . \
               matrix \
               scaling \
               fem \
               mesh \
#               extras \
               gui \
               material-data/can
#               material-data/freezing \
#               material-data/drying

# Input
HEADERS += arclength.h \
           basis.h \
           integrate.h \
           isoparam.h \
           mtxsolver.h \
           output.h \
#           extras/function.h \
#           extras/symbolic.h \
           scaling/scaling_ht.h \
           fem/auxsoln.h \
           fem/finite-element.h \
           fem/finite-element1d.h \
           fem/solution.h \
           gui/solver.h \
           gui/heat-gui.h \
           matrix/matrix.h \
           mesh/mesh1d.h \
           mesh/mesh2d.h \
           material-data/can/can.h \
           material-data/can/datafile.h
#           material-data/drying/drying.h \
#           material-data/drying/drying_D.h \
#           material-data/drying/regress.h \
#           material-data/freezing/freezing.h
FORMS += gui/solver.ui
SOURCES += basis.c \
           integrate.c \
           isoparam.c \
           mtxsolver.c \
           output.c \
           scaling/scaling_ht.c \
           fem/auxsoln.c \
           fem/finite-element.c \
           fem/finite-element1d.c \
           fem/solution.c \
           gui/main.cpp \
           gui/solver.cpp \
           gui/heat-gui.c \
           matrix/matrix.c \
           matrix/vector.c \
           mesh/mesh1d.c \
           mesh/mesh2d.c \
#           problems/2dlaplace.c \
#           problems/ce675p1.c \
#           problems/ce675p2.c \
#           problems/drying.c \
#           problems/fe-solver.c \
#           problems/heat-exact.c \
#           problems/heat-explicit.c \
#           problems/spheroid.c \
           material-data/can/can.c \
           material-data/can/datafile.c
#           material-data/drying/drying.c \
#           material-data/drying/drying_D.c \
#           material-data/drying/regress.c \
#           material-data/freezing/freezing.c
