#ifndef FREEZING_H
#define FREEZING_H

#include "ui_freezing.h"

extern "C" {
#include "matrix.h"
}

#include <qwt_plot_curve.h>

class Solver : public QMainWindow, private Ui::SolverWindow
{
    Q_OBJECT

    public:
        Solver(QWidget *parent = 0);
        ~Solver();
        
    public slots:
//        void solve();
        void quitApplication();
        void loadSimulation();
        void saveSimulation();
        void plotRandomProperty();

        void about();
        void aboutQt();
        void saveCSV();

        void changeLBC(int);
        void changeRBC(int);

        void setMaxNode(int);
        void setMaxTIndex(double);

        void solveProblems();
        void plotSolution();

    private:
        struct fe1d *problem;
        struct var *datalist;
        QwtPlotCurve *Temp, *Prod, *ice, *alpha, *prop;

        void leftBCHideAll();
        void rightBCHideAll();
        void leftBCShowConv();
        void leftBCShowConvTime();
        void rightBCShowConv();
        void rightBCShowConvTime();

        void storeVariables();
        void loadVars();

        void setupDomain();
        void plotResultsTime(int);
        void plotResultsSpace(int);
        void plotFunction(vector*, double (*)(double), QwtPlotCurve*);
};

#endif
