// ToDo: The composition data is now stored in a choi-okos structure instead of
// in a linked list; however, the code hasn't been tested yet. Make sure that
// it's not horribly broken.

#include <QtGui>

#include <qwt_plot.h> #include <qwt_plot_curve.h>
#include <qwt_math.h>

#include <stdio.h>
#include <string.h>

#include "solver.h"

extern "C" {
#include "heat-gui.h"
#include "solution.h"
#include "auxsoln.h"
#include "matrix.h"
#include "mesh1d.h"
#include "finite-element1d.h"
#include "solve.h"
#include "output.h"
#include "about.h"
}

#include "datafile.h"
#include "choi-okos.h"

extern double R;
choi_okos *comp_global;

Solver::Solver(QWidget *parent)
{
    setupUi(this); // Initialize the user interface

    // Connect all the relevant signals
    connect(buttonSolve, SIGNAL( clicked() ), this, SLOT( solveProblems() ));
    connect(buttonPlot, SIGNAL( clicked() ), this, SLOT( plotSolution() ));

    connect(actionQuit, SIGNAL( triggered() ), this, SLOT( quitApplication() ));
    connect(actionSave_Simulation, SIGNAL( triggered() ), this, SLOT( saveSimulation() ));
    connect(actionOpen, SIGNAL( triggered() ), this, SLOT( loadSimulation() ));
    connect(actionAbout, SIGNAL( triggered() ), this, SLOT( about() ));
    connect(actionAbout_Qt, SIGNAL( triggered() ), this, SLOT( aboutQt() ));
    connect(actionSave_CSV, SIGNAL( triggered() ), this, SLOT( saveCSV() ));

    connect(comboLeftBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeLBC(int) ));
    connect(comboRightBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeRBC(int) ));

    connect(spinNodes, SIGNAL( valueChanged(int) ), this, SLOT( setMaxNode(int) ));
    connect(spinDt, SIGNAL( valueChanged(double) ), this, SLOT(setMaxTIndex(double) ));
    connect(spintEnd, SIGNAL( valueChanged(double) ), this, SLOT(setMaxTIndex(double) ));

    // Set the domain and datalist variables to NULL. They're both initialized
    // later in the program.
    problem = NULL;
    datalist = NULL;

    // Create a structure to hold all the composition data
    composition = CreateChoiOkos();
    // Copy the address to a global variable so that it can be accessed from
    // the finite element solver.
    comp_global = composition;

    // Since the default boundary condition in the program is "insulation,"
    // none of the options need to be shown initially.
    leftBCHideAll();
    rightBCHideAll();

    // Check one of the radio buttons for the graph
    radioTime->setChecked(true);

    // Setup curves for each of the variables to plot. Also, make them
    // different colors.
    Temp = new QwtPlotCurve( "Temperature" );
    Temp->setPen( QPen( Qt::red ) );
    Prod = new QwtPlotCurve( "Product Concentration" );
    Prod->setPen( QPen( Qt::blue ) );
    Bact = new QwtPlotCurve( "Bacteria Concentration" );
    Bact->setPen( QPen( Qt::green ) );
    alpha = new QwtPlotCurve("Thermal Diffusivity");
    alpha->setPen( QPen( Qt::magenta ) );

    // Initialize the GUI
    comboLeftBC->setCurrentIndex(0);
    comboLeftBC->setEnabled(false);
    comboRightBC->setCurrentIndex(2);
    comboRightBC->setEnabled(false);
}

Solver::~Solver()
{
    DestroyFE1D(problem);
    delete Temp;
    delete Prod;
    delete Bact;
    delete alpha;

    DestroyChoiOkos(composition);
    comp_global = NULL;

    if(datalist)
        destroy_list(datalist);

    return;
}

// Close the window and exit the program.
void Solver::quitApplication()
{
    qApp->exit(0);
}

// Set the maximum value for the spin box that allows the user to select the
// node to plot data for.
void Solver::setMaxNode(int i) {
    spinPlotNode->setMaximum(i-1);
    return;
}

void Solver::setMaxTIndex(double i) {
    spinPlotTIndex->setMaximum((int) spintEnd->value()/spinDt->value()-1);
    return;
}

// When the user changes the boundary condition settings, show/hide options as
// appropriate.
void Solver::changeLBC(int index)
{
    // Hide everything to begin with.
    leftBCHideAll();
    switch(index) {
        case 1:
            // Show just options for convection.
            leftBCShowConv();
            break;
        case 2:
            // Show options for convection and heating/cooling.
            leftBCShowConvTime();
            break;
    }
    return;
}

// Same as above.
void Solver::changeRBC(int index)
{
    rightBCHideAll();
    switch(index) {
        case 1:
            rightBCShowConv();
            break;
        case 2:
            rightBCShowConvTime();
            break;
    }
    return;
}

// Function to hide all the widget associated with the left boundary condition.
void Solver::leftBCHideAll()
{
    lbch->hide();
    lbctext->hide();
    lbctexthot->hide();
    lbctextcold->hide();
    lbctheat->hide();
    spinHLeft->hide();
    spinTExtLeft->hide();
    spinTExtHotLeft->hide();
    spinTExtColdLeft->hide();
    spinTHeatLeft->hide();
    return;
}

// Same for the right BC.
void Solver::rightBCHideAll()
{
    rbch->hide();
    rbctext->hide();
    rbctexthot->hide();
    rbctextcold->hide();
    rbctheat->hide();
    spinHRight->hide();
    spinTExtRight->hide();
    spinTExtHotRight->hide();
    spinTExtColdRight->hide();
    spinTHeatRight->hide();

    return;
}

// Show options associated with convection.
void Solver::leftBCShowConv()
{
    lbch->show();
    lbctext->show();
    spinHLeft->show();
    spinTExtLeft->show();
 
    return;
}

void Solver::rightBCShowConv()
{
    rbch->show();
    rbctext->show();
    spinHRight->show();
    spinTExtRight->show();
    return;
}

// Show options for convection where the external temperature can take on two
// values.
void Solver::rightBCShowConvTime()
{
    rbch->show();
    rbctexthot->show();
    rbctextcold->show();
    rbctheat->show();
    spinHRight->show();
    spinTExtHotRight->show();
    spinTExtColdRight->show();
    spinTHeatRight->show();

    return;
}

void Solver::leftBCShowConvTime()
{
    lbch->show();
    lbctexthot->show();
    lbctextcold->show();
    lbctheat->show();
    spinHLeft->show();
    spinTExtHotLeft->show();
    spinTExtColdLeft->show();
    spinTHeatLeft->show();
    return;
}

// Macro to simplify the code for saving data from the gui to a linked list.
#define GETVARIABLE( VARNAME, BOXNAME ) { strcpy(tmp->name, #VARNAME); tmp->value = (BOXNAME)->value(); tmp->next = new_var(); tmp = tmp->next; }
// Store all (most of) the data that the user has inputted in the gui in the
// datalist variable.
void Solver::storeVariables()
{
    struct var *tmp;
    if(datalist)
        destroy_list(datalist);
    datalist = new_var();
    tmp = datalist;

    //GETVARIABLE(Mpro, spinMPro)
    //GETVARIABLE(Mfat, spinMFat)
    //GETVARIABLE(Mcar, spinMCar)
    //GETVARIABLE(Mfib, spinMFib)
    //GETVARIABLE(Mash, spinMAsh)
    //GETVARIABLE(Mwat, spinMWat)
    //GETVARIABLE(Mice, spinMIce)
    GETVARIABLE(AA, spinA1)
    GETVARIABLE(EaA, spinEa1)
    GETVARIABLE(AB, spinA2)
    GETVARIABLE(EaB, spinEa2)
    GETVARIABLE(To, spinTInit)
    GETVARIABLE(Text_hot, spinTExtHotRight)
    GETVARIABLE(Text_cold, spinTExtColdRight)
    GETVARIABLE(DomainWidth, spinWidth)
    GETVARIABLE(NNodes, spinNodes)
    GETVARIABLE(Deltat, spinDt)
    GETVARIABLE(t_heat, spinTHeatRight)
    GETVARIABLE(HConv, spinHRight)
    GETVARIABLE(CAinit, spinC1Init)
    GETVARIABLE(CBinit, spinC2Init)

    // Set the material composition here instead of in the linked list with the
    // rest of the values.
    composition->Mpro = spinMPro->value();
    composition->Mfat = spinMFat->value();
    composition->Mcar = spinMCar->value();
    composition->Mfib = spinMFib->value();
    composition->Mash = spinMAsh->value();
    composition->Mwat = spinMWat->value();
    composition->Mice = spinMIce->value();

    strcpy(tmp->name, "EndTime");
    tmp->value = spintEnd->value();

    return;
}

// Simplify code for moving values from a linked list into the gui.
#define RESTOREVAR( VARNAME, BOXNAME ) { if(strcmp(tmp->name, #VARNAME) == 0){(BOXNAME)->setValue(tmp->value);} }
// Take all the data from the datalist variable and stick it into the
// appropriate fields in the gui.
void Solver::loadVars()
{
    struct var *tmp;
    tmp = datalist;
    while(tmp) {
        //RESTOREVAR(Mpro, spinMPro)
        //RESTOREVAR(Mfat, spinMFat)
        //RESTOREVAR(Mcar, spinMCar)
        //RESTOREVAR(Mfib, spinMFib)
        //RESTOREVAR(Mash, spinMAsh)
        //RESTOREVAR(Mwat, spinMWat)
        //RESTOREVAR(Mice, spinMIce)
        RESTOREVAR(AA, spinA1)
        RESTOREVAR(EaA, spinEa1)
        RESTOREVAR(AB, spinA2)
        RESTOREVAR(EaB, spinEa2)
        RESTOREVAR(To, spinTInit)
        RESTOREVAR(Text_hot, spinTExtHotRight)
        RESTOREVAR(Text_cold, spinTExtColdRight)
        RESTOREVAR(NNodes, spinNodes)
        RESTOREVAR(Deltat, spinDt)
        RESTOREVAR(t_heat, spinTHeatRight)
        RESTOREVAR(HConv, spinHRight)
        RESTOREVAR(CAinit, spinC1Init)
        RESTOREVAR(CBinit, spinC2Init)
        RESTOREVAR(DomainWidth, spinWidth)
        RESTOREVAR(EndTime, spintEnd)

        if(strcmp(tmp->name, "R") == 0) {
            R = tmp->value;
        }

        tmp = tmp->next;
    }

    spinMPro->setValue(composition->Mpro);
    spinMFat->setValue(composition->Mfat);
    spinMCar->setValue(composition->Mcar);
    spinMFib->setValue(composition->Mfib);
    spinMAsh->setValue(composition->Mash);
    spinMWat->setValue(composition->Mwat);
    spinMIce->setValue(composition->Mice);

    return;
}

// Take the information about the domain from the appropriate fields in the gui
// and initialize the "domain" variable. Three dependant variables are solved
// for: Temperature, product concentration, and bacteria concentration.
void Solver::setupDomain()
{
    // Load the required parameters
    int NNodes = spinNodes->value();
    //double Deltax = spinWidth->value()/NNodes;
    double Deltat = spinDt->value();
    int NTimeSteps = (int) (spintEnd->value()/Deltat);

    Mesh1D *mesh;
    basis *b;
    matrix *IC;

    scaling_ht charvals; /* Scaling parameters */

    /* Define the characteristic temperature to be the initial temperature of
     * the can. */ 
    double Tc = spinTInit->value();
    /* Characteristic length = Domain radius */
    double Lc = spinWidth->value();
    /* Set the heat transfer coefficient as well. */
    double h = spinHRight->value();

    charvals = SetupScaling(alpha(Tc), Tc, Lc, k(Tc), h);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, 1.0, NNodes-1);

    /* Initialize most of the FE structure */
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         NTimeSteps);

    /* Set the values for dimensionless groups and stuff */
    problem->charvals = charvals;

    problem->nvars = 1; /* Set the number of variables to solve for */
    /* This bit might be a bit sketch */
    problem->dt = scaleTime(charvals, Deltat); /* Set the time step */

    /* Set the initial condition. Due to the way the dimensionless groups are
     * specified, the initial condition should always be 1. */
    IC = GenerateInitCondConst(problem, 0, 1);

    /* Initialize the transient solver */
    FE1DTransInit(problem, IC);

    /* Allocate space for solutions to two ODEs based on the PDE solution */
    FE1DInitAuxSolns(problem, 2);

    return;
}

void Solver::solveProblems()
{
    struct var *tmp;
    QProgressDialog progress("Solving PDE...", "Abort Solution", 0, 100, this);
    progress.setWindowModality(Qt::WindowModal);
    
    // Set all the global variables to 0 so that things don't break horribly
    // if everything wasn't defined.
    initialize_variables();

    // Set the datalist to NULL since we don't want to deal with old values.
    // This is a horrible idea, and the old data should be cleaned up instead.
    if(datalist)
        destroy_list(datalist);
    datalist = NULL;

    storeVariables();

    // Clean up stuff that was there before.
    DestroyFE1D(problem);

    // Save the data to the global variables used by the material property
    // funcitons. (Also bad.)
    tmp = datalist;
    while(tmp) {
        store_data(tmp);
        //if(strcmp(tmp->name, "HConv") == 0) {
        //    seth(tmp->value);
        //    printh();
        //}
        tmp = tmp->next;
    }

    // Initialize the problem
    setupDomain();
    progress.setMaximum(problem->maxsteps);

    print_global_vars();
    printf("rho = %g, k = %g, Cp = %g\n", rho(composition, 290), k(composition, 290), Cp(composition, 290));
    // Solve
    progress.setValue(0);
    while(problem->t<problem->maxsteps) {
        NLinSolve1DTransImp(problem, NULL);
        progress.setValue(problem->t);

        if(progress.wasCanceled()) {
            DestroyFE1D(problem);
            problem = NULL;
            break;
        }
    }
    
    if(problem) {
        progress.setValue(problem->maxsteps);

        SolveODE(problem, 0, 0, &react1, spinC1Init->value());
        SolveODE(problem, 0, 1, &react2, spinC2Init->value());

        plotSolution();
    }
}

void Solver::plotSolution()
{
    int nodenum = 0, timeindex = 0;

    // Check to see if the problem has already been solved. If not, report an
    // error and return before the program crashes.
    if(!problem) {
        QMessageBox::warning(this, "Error", "Please run the simulation before trying to plot stuff.");
        return;
    }

    if(!radioTime->isChecked()) {
        nodenum = spinPlotNode->value();
        plotResultsTime(nodenum);
    } else {
        timeindex = spinPlotTIndex->value();
        plotResultsSpace(timeindex);
    }
}

/* Example, since the Qwt docs are terrible

void Solver::plotResultsTime(struct Node1D* n)
{
    // Insert new curves
    QwtPlotCurve *cSin = new QwtPlotCurve( "y = sin(x)" );
    cSin->setPen( QPen( Qt::red ) ); // Draw it in red
    cSin->attach(qwtPlot); // Stick it in the plot

    // Create test data
    int npoints = 100;
    double *x, *y;
    x = (double*) calloc(npoints, sizeof(double));
    y = (double*) calloc(npoints, sizeof(double));

    int i;
    for(i=0; i<npoints; i++) {
        x[i] = i*3.14159/npoints;
        y[i] = sin(x[i]);
    }

    cSin->setRawSamples(x, y, npoints);

    // Show the plots
    qwtPlot->replot();
}
*/

// Function to plot results at a single node as a function of time.
void Solver::plotResultsTime(int nodenum)
{
    int i;
    int npts = problem->maxsteps;
    double *t, *T, *c, *d, *a, tmp;
    char *title;
    
    solution *s;

    title = (char*) calloc(20, sizeof(char));
    
    /* Create an array with the time values. */
    t = (double*) calloc(npts, sizeof(double));
    for(i=0; i<npts; i++) {
        t[i] = uscaleTime(problem->charvals, i*problem->dt);
    }

    // Detach the curves from the plot so that they don't show up. If we want
    // them to show up, we'll reattach them later.
    Temp->detach();
    Prod->detach();
    Bact->detach();
    alpha->detach();

    // Plot temperature data if the box for it is checked.
    if(checkTemp->isChecked()) {
        T = (double*) calloc(npts, sizeof(double));
        // Get the values from the solution.
        for(i=0; i<npts; i++) {
            s = FetchSolution(problem, i);
            T[i] = uscaleTemp(problem->charvals, val(s->val, nodenum, 0));
        }

        // Save the data to the curve. This makes a copy of the data, so we can
        // delete the two arrays later.
        Temp->setSamples(t, T, npts);

        // Put the curve on the graph.
        Temp->attach(qwtPlot);

        // Delete the temperature data we retrieved.
        free(T);
    }

    if(checkProduct->isChecked()) {
        c = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            s = FetchAuxSoln(problem, 0, i);
            c[i] = val(s->val, nodenum, 0);
        }

        Prod->attach(qwtPlot);

        Prod->setSamples(t, c, npts);
        free(c);
    }

    if(checkBact->isChecked()) {
        d = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            s = FetchAuxSoln(problem, 1, i);
            d[i] = val(s->val, nodenum, 0);
        }

        Bact->attach(qwtPlot);

        Bact->setSamples(t, d, npts);
        free(d);
    }
    
    if(checkk->isChecked()) {
        a = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            s = FetchSolution(problem, i);
            tmp = uscaleTemp(problem->charvals, val(s->val, nodenum, 0));
            a[i] = k(composition, tmp)/(rho(composition, tmp)*Cp(composition, tmp));
        }

        alpha->attach(qwtPlot);

        alpha->setSamples(t, a, npts);
        free(a);
    }

    sprintf(title, "x = %g",
            uscaleLength(problem->charvals,
                         valV(problem->mesh->nodes, nodenum)));
    qwtPlot->setTitle(title);
    free(title);

    // Set the axis title. 2 corresponds to the lower x axis.
    qwtPlot->setAxisTitle(2, "Time (sec)");

    // Update the graph to show what we just did.
    qwtPlot->replot();

    // Get rid of the "t" array. We don't want it anymore.
    free(t);
}

// Plot the results at a fixed time as a function of the spatial coordinate.
void Solver::plotResultsSpace(int t)
{
    int i;
    int npts = problem->nrows/problem->nvars;
    double *x, *T, *c, *d, *a, tmp;
    char *title;

    solution *s;

    title = (char*) calloc(20, sizeof(char));
    
    /* Create an array with the x values. */
    x = (double*) calloc(npts, sizeof(double));
    for(i=0; i<npts; i++) {
        /* This probably shouldn't be accessing the raw data from the vector */
        x[i] = uscaleLength(problem->charvals, valV(problem->mesh->nodes, i));
    }

    // Detach the curves from the plot so that they don't show up. If we want
    // them to show up, we'll reattach them later.
    Temp->detach();
    Prod->detach();
    Bact->detach();
    alpha->detach();

    // Plot temperature data if the box for it is checked.
    if(checkTemp->isChecked()) {
        T = (double*) calloc(npts, sizeof(double));
        // Get the values from the solution.
        s = FetchSolution(problem, t);
        for(i=0; i<npts; i++) {
            T[i] = uscaleTemp(problem->charvals, val(s->val, i, 0));
        }

        // Save the data to the curve. This makes a copy of the data, so we can
        // delete the two arrays later.
        Temp->setSamples(x, T, npts);

        // Put the curve on the graph.
        Temp->attach(qwtPlot);

        // Delete the temperature data we retrieved.
        free(T);
    }

    if(checkProduct->isChecked()) {
        c = (double*) calloc(npts, sizeof(double));
        s = FetchAuxSoln(problem, 0, t);
        for(i=0; i<npts; i++) {
            c[i] = val(s->val, i, 0);
        }

        Prod->attach(qwtPlot);

        Prod->setSamples(x, c, npts);
        free(c);
    }

    if(checkBact->isChecked()) {
        d = (double*) calloc(npts, sizeof(double));
        s = FetchAuxSoln(problem, 1, t);
        for(i=0; i<npts; i++) {
            d[i] = val(s->val, i, 0);
        }

        Bact->attach(qwtPlot);

        Bact->setSamples(x, d, npts);
        free(d);
    }
    
    if(checkk->isChecked()) {
        a = (double*) calloc(npts, sizeof(double));
        s = FetchSolution(problem, t);
        for(i=0; i<npts; i++) {
            tmp = uscaleTemp(problem->charvals, val(s->val, i, 0));

            a[i] = k(composition, tmp)/(rho(composition, tmp)*Cp(composition, tmp));
        }

        alpha->attach(qwtPlot);

        alpha->setSamples(x, a, npts);
        free(a);
    }

    // Set the title.
    sprintf(title, "t = %g", uscaleTime(problem->charvals, t*problem->dt));
    qwtPlot->setTitle(title);
    free(title);

    // Set the axis title. 2 corresponds to the lower x axis.
    qwtPlot->setAxisTitle(2, "Radius (m)");

    // Update the graph to show what we just did.
    qwtPlot->replot();

    // Don't get rid of the "x" array since it'll cause problems later.
    //free(x);
}

// Load all the simulation data from a text file.
void Solver::loadSimulation()
{
    struct var *tmp;
    QString path;
    QByteArray ba;
    char **buffer;
    int i;
    
    buffer = NULL;
    tmp = NULL;

    DestroyFE1D(problem);
    problem = NULL;
    if(datalist) {
        destroy_list(datalist);
        datalist = NULL;
    }

    // Get the filename to save to.
    path = QFileDialog::getOpenFileName(
            this,
            "Open",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();
    
    // Read the data file into a buffer.
    buffer = read_datafile(ba.data());

    // Parse the file line by line. (Hopefully there's not more than 200 lines.)
    for(i=0; i<200; i++) {
        // Strip any comments that might be lurking in the file.
        buffer[i] = remove_comments(buffer[i]);
        // Parse the line and save it to a variable.
        tmp = read_line(buffer[i]);
        // If the variable actually contained something, save it. Otherwise just
        // keep going and pretend like it never existed.
        // TODO: Fix the memory leak.
        if(strcmp(tmp->name, "NULL") != 0) {
            //printf("Name: %s -- Value: %g\n", tmp->name, tmp->value);
            datalist = push_var(datalist, tmp);
        }
    }

    // Clean up the buffer.
    delete_buffer(buffer);

    // Load the variables into the interface.
    loadVars();
    return;
}

// Save everything the user was working on to a file.
void Solver::saveSimulation()
{
    struct var *tmp;
    QString path;
    QByteArray ba;
    FILE *fp;
    storeVariables();
    tmp = datalist;


    // Get the filename to save to.
    path = QFileDialog::getSaveFileName(
            this,
            "Save As",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();

    // Open up the requested file.
    fp = fopen( (char*) ba.data(), "w+" );
    // Check for errors in opening the file.
    // TODO: Pop up an error message instead of printing stuff to the console.
    if(!fp) {
        fprintf(stderr, "Failed to open file for writing.");
        return;
    }

    while(tmp) {
        fprintf(fp, "%s = %g\n", tmp->name, tmp->value);
        tmp = tmp->next;
    }

    // Close the file.
    fclose(fp);

    return;
}

void Solver::saveCSV()
{
    QString path;
    QByteArray ba;

    // Check to see if the user has run the simulation before trying to save
    // the output. Stops the program from seg faulting.
    if(!problem) {
        QMessageBox::warning(this, "Error", "Please run the simulation before saving the results.");
        return;
    }

    // Get the filename to save to.
    path = QFileDialog::getSaveFileName(
            this,
            "Save As",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();

    if(radioTime->isChecked()) {
        CSVOutFixedTime(problem, spinPlotTIndex->value(), ba.data());
    } else {
        CSVOutFixedNode(problem, spinPlotNode->value(), ba.data());
    }
    return;
}

// Bring up an about dialog.
void Solver::about()
{
    QMessageBox::about(this, "About", genAboutString());
}

void Solver::aboutQt()
{
    QMessageBox::aboutQt(this, "About Qt");
}

