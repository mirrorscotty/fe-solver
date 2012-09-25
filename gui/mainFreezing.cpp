#include <QApplication>

#include "solver.h"

// Main loop for the gui.
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Solver *main_window = new Solver;

    main_window->show();

    return app.exec();
}

