#ifndef SOLVER_H
#define SOLVER_H

#include "ui_solver.h"
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
        choi_okos *composition;
        QwtPlotCurve *Temp, *Prod, *Bact, *alpha;

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
};

#endif
