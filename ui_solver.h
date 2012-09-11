/********************************************************************************
** Form generated from reading UI file 'solver.ui'
**
** Created: Tue Sep 11 02:10:48 2012
**      by: Qt User Interface Compiler version 4.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SOLVER_H
#define UI_SOLVER_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_SolverWindow
{
public:
    QAction *actionQuit;
    QAction *actionAbout;
    QAction *actionSave_Simulation;
    QAction *actionOpen;
    QAction *actionSave_CSV;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout_29;
    QToolBox *toolBox;
    QWidget *DomainTab;
    QVBoxLayout *verticalLayout_8;
    QLabel *label_37;
    QFrame *line_11;
    QHBoxLayout *horizontalLayout_33;
    QLabel *label_38;
    QSpinBox *spinNodes;
    QHBoxLayout *horizontalLayout_32;
    QLabel *label_39;
    QDoubleSpinBox *spinWidth;
    QLabel *label_40;
    QFrame *line_12;
    QHBoxLayout *horizontalLayout_31;
    QLabel *label_41;
    QDoubleSpinBox *spinDt;
    QHBoxLayout *horizontalLayout_30;
    QLabel *label_42;
    QDoubleSpinBox *spintEnd;
    QSpacerItem *verticalSpacer_7;
    QWidget *MatTab;
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QFrame *line;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_2;
    QDoubleSpinBox *spinMPro;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_3;
    QDoubleSpinBox *spinMFat;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_4;
    QDoubleSpinBox *spinMCar;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_5;
    QDoubleSpinBox *spinMFib;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_6;
    QDoubleSpinBox *spinMAsh;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_7;
    QDoubleSpinBox *spinMWat;
    QHBoxLayout *horizontalLayout;
    QLabel *label_8;
    QDoubleSpinBox *spinMIce;
    QSpacerItem *verticalSpacer;
    QWidget *TTab;
    QVBoxLayout *verticalLayout_4;
    QLabel *label_11;
    QFrame *line_4;
    QHBoxLayout *horizontalLayout_10;
    QLabel *label_12;
    QDoubleSpinBox *spinTInit;
    QLabel *label_13;
    QFrame *line_5;
    QHBoxLayout *horizontalLayout_15;
    QLabel *label_15;
    QComboBox *comboLeftBC;
    QHBoxLayout *horizontalLayout_11;
    QLabel *lbch;
    QDoubleSpinBox *spinHLeft;
    QHBoxLayout *horizontalLayout_14;
    QLabel *lbctext;
    QDoubleSpinBox *spinTExtLeft;
    QHBoxLayout *horizontalLayout_13;
    QLabel *lbctexthot;
    QDoubleSpinBox *spinTExtHotLeft;
    QHBoxLayout *horizontalLayout_12;
    QLabel *lbctextcold;
    QDoubleSpinBox *spinTExtColdLeft;
    QHBoxLayout *horizontalLayout_20;
    QLabel *lbctheat;
    QDoubleSpinBox *spinTHeatLeft;
    QLabel *label_19;
    QFrame *line_6;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_14;
    QComboBox *comboRightBC;
    QHBoxLayout *horizontalLayout_16;
    QLabel *rbch;
    QDoubleSpinBox *spinHRight;
    QHBoxLayout *horizontalLayout_17;
    QLabel *rbctext;
    QDoubleSpinBox *spinTExtRight;
    QHBoxLayout *horizontalLayout_18;
    QLabel *rbctexthot;
    QDoubleSpinBox *spinTExtHotRight;
    QHBoxLayout *horizontalLayout_19;
    QLabel *rbctextcold;
    QDoubleSpinBox *spinTExtColdRight;
    QHBoxLayout *horizontalLayout_21;
    QLabel *rbctheat;
    QDoubleSpinBox *spinTHeatRight;
    QSpacerItem *verticalSpacer_2;
    QWidget *C1Tab;
    QVBoxLayout *verticalLayout_5;
    QLabel *label_31;
    QFrame *line_8;
    QHBoxLayout *horizontalLayout_24;
    QLabel *label_30;
    QDoubleSpinBox *spinC1Init;
    QLabel *label_28;
    QFrame *line_7;
    QHBoxLayout *horizontalLayout_23;
    QLabel *label_27;
    QDoubleSpinBox *spinA1;
    QHBoxLayout *horizontalLayout_22;
    QLabel *label_29;
    QDoubleSpinBox *spinEa1;
    QSpacerItem *verticalSpacer_3;
    QWidget *C2Tab;
    QVBoxLayout *verticalLayout_6;
    QLabel *label_33;
    QFrame *line_10;
    QHBoxLayout *horizontalLayout_27;
    QLabel *label_35;
    QDoubleSpinBox *spinC2Init;
    QLabel *label_32;
    QFrame *line_9;
    QHBoxLayout *horizontalLayout_26;
    QLabel *label_34;
    QDoubleSpinBox *spinA2;
    QHBoxLayout *horizontalLayout_25;
    QLabel *label_36;
    QDoubleSpinBox *spinEa2;
    QSpacerItem *verticalSpacer_4;
    QVBoxLayout *verticalLayout_7;
    QwtPlot *qwtPlot;
    QHBoxLayout *horizontalLayout_28;
    QGroupBox *DepGroup;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *checkTemp;
    QCheckBox *checkk;
    QCheckBox *checkProduct;
    QCheckBox *checkBact;
    QGroupBox *IndepGroup;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_34;
    QRadioButton *radioTime;
    QLabel *label_9;
    QSpinBox *spinPlotTIndex;
    QHBoxLayout *horizontalLayout_35;
    QRadioButton *radioLength;
    QLabel *label_10;
    QSpinBox *spinPlotNode;
    QPushButton *buttonPlot;
    QHBoxLayout *horizontalLayout_8;
    QProgressBar *progressBar;
    QPushButton *buttonSolve;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *SolverWindow)
    {
        if (SolverWindow->objectName().isEmpty())
            SolverWindow->setObjectName(QString::fromUtf8("SolverWindow"));
        SolverWindow->resize(826, 680);
        SolverWindow->setAutoFillBackground(false);
        actionQuit = new QAction(SolverWindow);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionAbout = new QAction(SolverWindow);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionSave_Simulation = new QAction(SolverWindow);
        actionSave_Simulation->setObjectName(QString::fromUtf8("actionSave_Simulation"));
        actionOpen = new QAction(SolverWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionSave_CSV = new QAction(SolverWindow);
        actionSave_CSV->setObjectName(QString::fromUtf8("actionSave_CSV"));
        centralwidget = new QWidget(SolverWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        horizontalLayout_29 = new QHBoxLayout(centralwidget);
        horizontalLayout_29->setObjectName(QString::fromUtf8("horizontalLayout_29"));
        toolBox = new QToolBox(centralwidget);
        toolBox->setObjectName(QString::fromUtf8("toolBox"));
        toolBox->setEnabled(true);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(toolBox->sizePolicy().hasHeightForWidth());
        toolBox->setSizePolicy(sizePolicy);
        toolBox->setMinimumSize(QSize(350, 0));
        toolBox->setMaximumSize(QSize(350, 16777215));
        DomainTab = new QWidget();
        DomainTab->setObjectName(QString::fromUtf8("DomainTab"));
        DomainTab->setGeometry(QRect(0, 0, 340, 493));
        verticalLayout_8 = new QVBoxLayout(DomainTab);
        verticalLayout_8->setObjectName(QString::fromUtf8("verticalLayout_8"));
        label_37 = new QLabel(DomainTab);
        label_37->setObjectName(QString::fromUtf8("label_37"));

        verticalLayout_8->addWidget(label_37);

        line_11 = new QFrame(DomainTab);
        line_11->setObjectName(QString::fromUtf8("line_11"));
        line_11->setFrameShape(QFrame::HLine);
        line_11->setFrameShadow(QFrame::Sunken);

        verticalLayout_8->addWidget(line_11);

        horizontalLayout_33 = new QHBoxLayout();
        horizontalLayout_33->setObjectName(QString::fromUtf8("horizontalLayout_33"));
        label_38 = new QLabel(DomainTab);
        label_38->setObjectName(QString::fromUtf8("label_38"));

        horizontalLayout_33->addWidget(label_38);

        spinNodes = new QSpinBox(DomainTab);
        spinNodes->setObjectName(QString::fromUtf8("spinNodes"));
        spinNodes->setMinimum(2);
        spinNodes->setMaximum(999);

        horizontalLayout_33->addWidget(spinNodes);


        verticalLayout_8->addLayout(horizontalLayout_33);

        horizontalLayout_32 = new QHBoxLayout();
        horizontalLayout_32->setObjectName(QString::fromUtf8("horizontalLayout_32"));
        label_39 = new QLabel(DomainTab);
        label_39->setObjectName(QString::fromUtf8("label_39"));

        horizontalLayout_32->addWidget(label_39);

        spinWidth = new QDoubleSpinBox(DomainTab);
        spinWidth->setObjectName(QString::fromUtf8("spinWidth"));
        spinWidth->setMinimum(0.1);

        horizontalLayout_32->addWidget(spinWidth);


        verticalLayout_8->addLayout(horizontalLayout_32);

        label_40 = new QLabel(DomainTab);
        label_40->setObjectName(QString::fromUtf8("label_40"));

        verticalLayout_8->addWidget(label_40);

        line_12 = new QFrame(DomainTab);
        line_12->setObjectName(QString::fromUtf8("line_12"));
        line_12->setFrameShape(QFrame::HLine);
        line_12->setFrameShadow(QFrame::Sunken);

        verticalLayout_8->addWidget(line_12);

        horizontalLayout_31 = new QHBoxLayout();
        horizontalLayout_31->setObjectName(QString::fromUtf8("horizontalLayout_31"));
        label_41 = new QLabel(DomainTab);
        label_41->setObjectName(QString::fromUtf8("label_41"));

        horizontalLayout_31->addWidget(label_41);

        spinDt = new QDoubleSpinBox(DomainTab);
        spinDt->setObjectName(QString::fromUtf8("spinDt"));
        spinDt->setDecimals(4);
        spinDt->setMinimum(0.0001);

        horizontalLayout_31->addWidget(spinDt);


        verticalLayout_8->addLayout(horizontalLayout_31);

        horizontalLayout_30 = new QHBoxLayout();
        horizontalLayout_30->setObjectName(QString::fromUtf8("horizontalLayout_30"));
        label_42 = new QLabel(DomainTab);
        label_42->setObjectName(QString::fromUtf8("label_42"));

        horizontalLayout_30->addWidget(label_42);

        spintEnd = new QDoubleSpinBox(DomainTab);
        spintEnd->setObjectName(QString::fromUtf8("spintEnd"));
        spintEnd->setMinimum(0.1);
        spintEnd->setMaximum(1e+06);

        horizontalLayout_30->addWidget(spintEnd);


        verticalLayout_8->addLayout(horizontalLayout_30);

        verticalSpacer_7 = new QSpacerItem(20, 284, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_8->addItem(verticalSpacer_7);

        toolBox->addItem(DomainTab, QString::fromUtf8("Domain Settings"));
        MatTab = new QWidget();
        MatTab->setObjectName(QString::fromUtf8("MatTab"));
        MatTab->setGeometry(QRect(0, 0, 340, 493));
        verticalLayout = new QVBoxLayout(MatTab);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        label = new QLabel(MatTab);
        label->setObjectName(QString::fromUtf8("label"));
        QFont font;
        font.setPointSize(12);
        label->setFont(font);

        verticalLayout->addWidget(label);

        line = new QFrame(MatTab);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        label_2 = new QLabel(MatTab);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_7->addWidget(label_2);

        spinMPro = new QDoubleSpinBox(MatTab);
        spinMPro->setObjectName(QString::fromUtf8("spinMPro"));
        spinMPro->setMaximum(1);
        spinMPro->setSingleStep(0.1);

        horizontalLayout_7->addWidget(spinMPro);


        verticalLayout->addLayout(horizontalLayout_7);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        label_3 = new QLabel(MatTab);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_6->addWidget(label_3);

        spinMFat = new QDoubleSpinBox(MatTab);
        spinMFat->setObjectName(QString::fromUtf8("spinMFat"));
        spinMFat->setMaximum(1);
        spinMFat->setSingleStep(0.1);

        horizontalLayout_6->addWidget(spinMFat);


        verticalLayout->addLayout(horizontalLayout_6);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        label_4 = new QLabel(MatTab);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        horizontalLayout_5->addWidget(label_4);

        spinMCar = new QDoubleSpinBox(MatTab);
        spinMCar->setObjectName(QString::fromUtf8("spinMCar"));
        spinMCar->setMaximum(1);
        spinMCar->setSingleStep(0.1);

        horizontalLayout_5->addWidget(spinMCar);


        verticalLayout->addLayout(horizontalLayout_5);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        label_5 = new QLabel(MatTab);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        horizontalLayout_4->addWidget(label_5);

        spinMFib = new QDoubleSpinBox(MatTab);
        spinMFib->setObjectName(QString::fromUtf8("spinMFib"));
        spinMFib->setMaximum(1);
        spinMFib->setSingleStep(0.1);

        horizontalLayout_4->addWidget(spinMFib);


        verticalLayout->addLayout(horizontalLayout_4);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label_6 = new QLabel(MatTab);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        horizontalLayout_3->addWidget(label_6);

        spinMAsh = new QDoubleSpinBox(MatTab);
        spinMAsh->setObjectName(QString::fromUtf8("spinMAsh"));
        spinMAsh->setMaximum(1);
        spinMAsh->setSingleStep(0.1);

        horizontalLayout_3->addWidget(spinMAsh);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_7 = new QLabel(MatTab);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        horizontalLayout_2->addWidget(label_7);

        spinMWat = new QDoubleSpinBox(MatTab);
        spinMWat->setObjectName(QString::fromUtf8("spinMWat"));
        spinMWat->setMaximum(1);
        spinMWat->setSingleStep(0.1);

        horizontalLayout_2->addWidget(spinMWat);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label_8 = new QLabel(MatTab);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        horizontalLayout->addWidget(label_8);

        spinMIce = new QDoubleSpinBox(MatTab);
        spinMIce->setObjectName(QString::fromUtf8("spinMIce"));
        spinMIce->setMaximum(1);
        spinMIce->setSingleStep(0.1);

        horizontalLayout->addWidget(spinMIce);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer = new QSpacerItem(20, 254, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        toolBox->addItem(MatTab, QString::fromUtf8("Material Properties"));
        line->raise();
        label->raise();
        TTab = new QWidget();
        TTab->setObjectName(QString::fromUtf8("TTab"));
        TTab->setGeometry(QRect(0, 0, 340, 493));
        verticalLayout_4 = new QVBoxLayout(TTab);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        label_11 = new QLabel(TTab);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        verticalLayout_4->addWidget(label_11);

        line_4 = new QFrame(TTab);
        line_4->setObjectName(QString::fromUtf8("line_4"));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line_4);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        label_12 = new QLabel(TTab);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        horizontalLayout_10->addWidget(label_12);

        spinTInit = new QDoubleSpinBox(TTab);
        spinTInit->setObjectName(QString::fromUtf8("spinTInit"));
        spinTInit->setMaximum(9999.99);

        horizontalLayout_10->addWidget(spinTInit);


        verticalLayout_4->addLayout(horizontalLayout_10);

        label_13 = new QLabel(TTab);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        verticalLayout_4->addWidget(label_13);

        line_5 = new QFrame(TTab);
        line_5->setObjectName(QString::fromUtf8("line_5"));
        line_5->setFrameShape(QFrame::HLine);
        line_5->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line_5);

        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setObjectName(QString::fromUtf8("horizontalLayout_15"));
        label_15 = new QLabel(TTab);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        horizontalLayout_15->addWidget(label_15);

        comboLeftBC = new QComboBox(TTab);
        comboLeftBC->setObjectName(QString::fromUtf8("comboLeftBC"));

        horizontalLayout_15->addWidget(comboLeftBC);


        verticalLayout_4->addLayout(horizontalLayout_15);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        lbch = new QLabel(TTab);
        lbch->setObjectName(QString::fromUtf8("lbch"));

        horizontalLayout_11->addWidget(lbch);

        spinHLeft = new QDoubleSpinBox(TTab);
        spinHLeft->setObjectName(QString::fromUtf8("spinHLeft"));
        spinHLeft->setMaximum(999.99);

        horizontalLayout_11->addWidget(spinHLeft);


        verticalLayout_4->addLayout(horizontalLayout_11);

        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setObjectName(QString::fromUtf8("horizontalLayout_14"));
        lbctext = new QLabel(TTab);
        lbctext->setObjectName(QString::fromUtf8("lbctext"));

        horizontalLayout_14->addWidget(lbctext);

        spinTExtLeft = new QDoubleSpinBox(TTab);
        spinTExtLeft->setObjectName(QString::fromUtf8("spinTExtLeft"));
        spinTExtLeft->setMaximum(9999.99);

        horizontalLayout_14->addWidget(spinTExtLeft);


        verticalLayout_4->addLayout(horizontalLayout_14);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        lbctexthot = new QLabel(TTab);
        lbctexthot->setObjectName(QString::fromUtf8("lbctexthot"));

        horizontalLayout_13->addWidget(lbctexthot);

        spinTExtHotLeft = new QDoubleSpinBox(TTab);
        spinTExtHotLeft->setObjectName(QString::fromUtf8("spinTExtHotLeft"));
        spinTExtHotLeft->setMaximum(9999.99);

        horizontalLayout_13->addWidget(spinTExtHotLeft);


        verticalLayout_4->addLayout(horizontalLayout_13);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        lbctextcold = new QLabel(TTab);
        lbctextcold->setObjectName(QString::fromUtf8("lbctextcold"));

        horizontalLayout_12->addWidget(lbctextcold);

        spinTExtColdLeft = new QDoubleSpinBox(TTab);
        spinTExtColdLeft->setObjectName(QString::fromUtf8("spinTExtColdLeft"));
        spinTExtColdLeft->setMaximum(9999.99);

        horizontalLayout_12->addWidget(spinTExtColdLeft);


        verticalLayout_4->addLayout(horizontalLayout_12);

        horizontalLayout_20 = new QHBoxLayout();
        horizontalLayout_20->setObjectName(QString::fromUtf8("horizontalLayout_20"));
        lbctheat = new QLabel(TTab);
        lbctheat->setObjectName(QString::fromUtf8("lbctheat"));

        horizontalLayout_20->addWidget(lbctheat);

        spinTHeatLeft = new QDoubleSpinBox(TTab);
        spinTHeatLeft->setObjectName(QString::fromUtf8("spinTHeatLeft"));
        spinTHeatLeft->setMaximum(1e+06);

        horizontalLayout_20->addWidget(spinTHeatLeft);


        verticalLayout_4->addLayout(horizontalLayout_20);

        label_19 = new QLabel(TTab);
        label_19->setObjectName(QString::fromUtf8("label_19"));

        verticalLayout_4->addWidget(label_19);

        line_6 = new QFrame(TTab);
        line_6->setObjectName(QString::fromUtf8("line_6"));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);

        verticalLayout_4->addWidget(line_6);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        label_14 = new QLabel(TTab);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        horizontalLayout_9->addWidget(label_14);

        comboRightBC = new QComboBox(TTab);
        comboRightBC->setObjectName(QString::fromUtf8("comboRightBC"));

        horizontalLayout_9->addWidget(comboRightBC);


        verticalLayout_4->addLayout(horizontalLayout_9);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setObjectName(QString::fromUtf8("horizontalLayout_16"));
        rbch = new QLabel(TTab);
        rbch->setObjectName(QString::fromUtf8("rbch"));

        horizontalLayout_16->addWidget(rbch);

        spinHRight = new QDoubleSpinBox(TTab);
        spinHRight->setObjectName(QString::fromUtf8("spinHRight"));
        spinHRight->setMaximum(9999.99);

        horizontalLayout_16->addWidget(spinHRight);


        verticalLayout_4->addLayout(horizontalLayout_16);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        rbctext = new QLabel(TTab);
        rbctext->setObjectName(QString::fromUtf8("rbctext"));

        horizontalLayout_17->addWidget(rbctext);

        spinTExtRight = new QDoubleSpinBox(TTab);
        spinTExtRight->setObjectName(QString::fromUtf8("spinTExtRight"));
        spinTExtRight->setMaximum(9999.99);

        horizontalLayout_17->addWidget(spinTExtRight);


        verticalLayout_4->addLayout(horizontalLayout_17);

        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setObjectName(QString::fromUtf8("horizontalLayout_18"));
        rbctexthot = new QLabel(TTab);
        rbctexthot->setObjectName(QString::fromUtf8("rbctexthot"));

        horizontalLayout_18->addWidget(rbctexthot);

        spinTExtHotRight = new QDoubleSpinBox(TTab);
        spinTExtHotRight->setObjectName(QString::fromUtf8("spinTExtHotRight"));
        spinTExtHotRight->setMaximum(9999.99);

        horizontalLayout_18->addWidget(spinTExtHotRight);


        verticalLayout_4->addLayout(horizontalLayout_18);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setObjectName(QString::fromUtf8("horizontalLayout_19"));
        rbctextcold = new QLabel(TTab);
        rbctextcold->setObjectName(QString::fromUtf8("rbctextcold"));

        horizontalLayout_19->addWidget(rbctextcold);

        spinTExtColdRight = new QDoubleSpinBox(TTab);
        spinTExtColdRight->setObjectName(QString::fromUtf8("spinTExtColdRight"));
        spinTExtColdRight->setMaximum(9999.99);

        horizontalLayout_19->addWidget(spinTExtColdRight);


        verticalLayout_4->addLayout(horizontalLayout_19);

        horizontalLayout_21 = new QHBoxLayout();
        horizontalLayout_21->setObjectName(QString::fromUtf8("horizontalLayout_21"));
        rbctheat = new QLabel(TTab);
        rbctheat->setObjectName(QString::fromUtf8("rbctheat"));

        horizontalLayout_21->addWidget(rbctheat);

        spinTHeatRight = new QDoubleSpinBox(TTab);
        spinTHeatRight->setObjectName(QString::fromUtf8("spinTHeatRight"));
        spinTHeatRight->setMaximum(1e+06);

        horizontalLayout_21->addWidget(spinTHeatRight);


        verticalLayout_4->addLayout(horizontalLayout_21);

        verticalSpacer_2 = new QSpacerItem(20, 8, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer_2);

        toolBox->addItem(TTab, QString::fromUtf8("Heat Conduction"));
        C1Tab = new QWidget();
        C1Tab->setObjectName(QString::fromUtf8("C1Tab"));
        C1Tab->setGeometry(QRect(0, 0, 340, 493));
        verticalLayout_5 = new QVBoxLayout(C1Tab);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        label_31 = new QLabel(C1Tab);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        verticalLayout_5->addWidget(label_31);

        line_8 = new QFrame(C1Tab);
        line_8->setObjectName(QString::fromUtf8("line_8"));
        line_8->setFrameShape(QFrame::HLine);
        line_8->setFrameShadow(QFrame::Sunken);

        verticalLayout_5->addWidget(line_8);

        horizontalLayout_24 = new QHBoxLayout();
        horizontalLayout_24->setObjectName(QString::fromUtf8("horizontalLayout_24"));
        label_30 = new QLabel(C1Tab);
        label_30->setObjectName(QString::fromUtf8("label_30"));

        horizontalLayout_24->addWidget(label_30);

        spinC1Init = new QDoubleSpinBox(C1Tab);
        spinC1Init->setObjectName(QString::fromUtf8("spinC1Init"));

        horizontalLayout_24->addWidget(spinC1Init);


        verticalLayout_5->addLayout(horizontalLayout_24);

        label_28 = new QLabel(C1Tab);
        label_28->setObjectName(QString::fromUtf8("label_28"));

        verticalLayout_5->addWidget(label_28);

        line_7 = new QFrame(C1Tab);
        line_7->setObjectName(QString::fromUtf8("line_7"));
        line_7->setFrameShape(QFrame::HLine);
        line_7->setFrameShadow(QFrame::Sunken);

        verticalLayout_5->addWidget(line_7);

        horizontalLayout_23 = new QHBoxLayout();
        horizontalLayout_23->setObjectName(QString::fromUtf8("horizontalLayout_23"));
        label_27 = new QLabel(C1Tab);
        label_27->setObjectName(QString::fromUtf8("label_27"));

        horizontalLayout_23->addWidget(label_27);

        spinA1 = new QDoubleSpinBox(C1Tab);
        spinA1->setObjectName(QString::fromUtf8("spinA1"));

        horizontalLayout_23->addWidget(spinA1);


        verticalLayout_5->addLayout(horizontalLayout_23);

        horizontalLayout_22 = new QHBoxLayout();
        horizontalLayout_22->setObjectName(QString::fromUtf8("horizontalLayout_22"));
        label_29 = new QLabel(C1Tab);
        label_29->setObjectName(QString::fromUtf8("label_29"));

        horizontalLayout_22->addWidget(label_29);

        spinEa1 = new QDoubleSpinBox(C1Tab);
        spinEa1->setObjectName(QString::fromUtf8("spinEa1"));
        spinEa1->setMaximum(1e+07);

        horizontalLayout_22->addWidget(spinEa1);


        verticalLayout_5->addLayout(horizontalLayout_22);

        verticalSpacer_3 = new QSpacerItem(20, 342, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_5->addItem(verticalSpacer_3);

        toolBox->addItem(C1Tab, QString::fromUtf8("Product Reaction Rate"));
        C2Tab = new QWidget();
        C2Tab->setObjectName(QString::fromUtf8("C2Tab"));
        C2Tab->setGeometry(QRect(0, 0, 340, 493));
        verticalLayout_6 = new QVBoxLayout(C2Tab);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        label_33 = new QLabel(C2Tab);
        label_33->setObjectName(QString::fromUtf8("label_33"));

        verticalLayout_6->addWidget(label_33);

        line_10 = new QFrame(C2Tab);
        line_10->setObjectName(QString::fromUtf8("line_10"));
        line_10->setFrameShape(QFrame::HLine);
        line_10->setFrameShadow(QFrame::Sunken);

        verticalLayout_6->addWidget(line_10);

        horizontalLayout_27 = new QHBoxLayout();
        horizontalLayout_27->setObjectName(QString::fromUtf8("horizontalLayout_27"));
        label_35 = new QLabel(C2Tab);
        label_35->setObjectName(QString::fromUtf8("label_35"));

        horizontalLayout_27->addWidget(label_35);

        spinC2Init = new QDoubleSpinBox(C2Tab);
        spinC2Init->setObjectName(QString::fromUtf8("spinC2Init"));

        horizontalLayout_27->addWidget(spinC2Init);


        verticalLayout_6->addLayout(horizontalLayout_27);

        label_32 = new QLabel(C2Tab);
        label_32->setObjectName(QString::fromUtf8("label_32"));

        verticalLayout_6->addWidget(label_32);

        line_9 = new QFrame(C2Tab);
        line_9->setObjectName(QString::fromUtf8("line_9"));
        line_9->setFrameShape(QFrame::HLine);
        line_9->setFrameShadow(QFrame::Sunken);

        verticalLayout_6->addWidget(line_9);

        horizontalLayout_26 = new QHBoxLayout();
        horizontalLayout_26->setObjectName(QString::fromUtf8("horizontalLayout_26"));
        label_34 = new QLabel(C2Tab);
        label_34->setObjectName(QString::fromUtf8("label_34"));

        horizontalLayout_26->addWidget(label_34);

        spinA2 = new QDoubleSpinBox(C2Tab);
        spinA2->setObjectName(QString::fromUtf8("spinA2"));

        horizontalLayout_26->addWidget(spinA2);


        verticalLayout_6->addLayout(horizontalLayout_26);

        horizontalLayout_25 = new QHBoxLayout();
        horizontalLayout_25->setObjectName(QString::fromUtf8("horizontalLayout_25"));
        label_36 = new QLabel(C2Tab);
        label_36->setObjectName(QString::fromUtf8("label_36"));

        horizontalLayout_25->addWidget(label_36);

        spinEa2 = new QDoubleSpinBox(C2Tab);
        spinEa2->setObjectName(QString::fromUtf8("spinEa2"));
        spinEa2->setMaximum(1e+07);

        horizontalLayout_25->addWidget(spinEa2);


        verticalLayout_6->addLayout(horizontalLayout_25);

        verticalSpacer_4 = new QSpacerItem(20, 342, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_6->addItem(verticalSpacer_4);

        toolBox->addItem(C2Tab, QString::fromUtf8("Bacteria Reaction Rate"));

        horizontalLayout_29->addWidget(toolBox);

        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        qwtPlot = new QwtPlot(centralwidget);
        qwtPlot->setObjectName(QString::fromUtf8("qwtPlot"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(qwtPlot->sizePolicy().hasHeightForWidth());
        qwtPlot->setSizePolicy(sizePolicy1);

        verticalLayout_7->addWidget(qwtPlot);

        horizontalLayout_28 = new QHBoxLayout();
        horizontalLayout_28->setObjectName(QString::fromUtf8("horizontalLayout_28"));
        DepGroup = new QGroupBox(centralwidget);
        DepGroup->setObjectName(QString::fromUtf8("DepGroup"));
        verticalLayout_3 = new QVBoxLayout(DepGroup);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        checkTemp = new QCheckBox(DepGroup);
        checkTemp->setObjectName(QString::fromUtf8("checkTemp"));

        verticalLayout_3->addWidget(checkTemp);

        checkk = new QCheckBox(DepGroup);
        checkk->setObjectName(QString::fromUtf8("checkk"));
        checkk->setLayoutDirection(Qt::LeftToRight);

        verticalLayout_3->addWidget(checkk);

        checkProduct = new QCheckBox(DepGroup);
        checkProduct->setObjectName(QString::fromUtf8("checkProduct"));

        verticalLayout_3->addWidget(checkProduct);

        checkBact = new QCheckBox(DepGroup);
        checkBact->setObjectName(QString::fromUtf8("checkBact"));

        verticalLayout_3->addWidget(checkBact);


        horizontalLayout_28->addWidget(DepGroup);

        IndepGroup = new QGroupBox(centralwidget);
        IndepGroup->setObjectName(QString::fromUtf8("IndepGroup"));
        verticalLayout_2 = new QVBoxLayout(IndepGroup);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout_34 = new QHBoxLayout();
        horizontalLayout_34->setObjectName(QString::fromUtf8("horizontalLayout_34"));
        radioTime = new QRadioButton(IndepGroup);
        radioTime->setObjectName(QString::fromUtf8("radioTime"));

        horizontalLayout_34->addWidget(radioTime);

        label_9 = new QLabel(IndepGroup);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        horizontalLayout_34->addWidget(label_9);

        spinPlotTIndex = new QSpinBox(IndepGroup);
        spinPlotTIndex->setObjectName(QString::fromUtf8("spinPlotTIndex"));

        horizontalLayout_34->addWidget(spinPlotTIndex);


        verticalLayout_2->addLayout(horizontalLayout_34);

        horizontalLayout_35 = new QHBoxLayout();
        horizontalLayout_35->setObjectName(QString::fromUtf8("horizontalLayout_35"));
        radioLength = new QRadioButton(IndepGroup);
        radioLength->setObjectName(QString::fromUtf8("radioLength"));

        horizontalLayout_35->addWidget(radioLength);

        label_10 = new QLabel(IndepGroup);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        horizontalLayout_35->addWidget(label_10);

        spinPlotNode = new QSpinBox(IndepGroup);
        spinPlotNode->setObjectName(QString::fromUtf8("spinPlotNode"));

        horizontalLayout_35->addWidget(spinPlotNode);


        verticalLayout_2->addLayout(horizontalLayout_35);

        buttonPlot = new QPushButton(IndepGroup);
        buttonPlot->setObjectName(QString::fromUtf8("buttonPlot"));

        verticalLayout_2->addWidget(buttonPlot);


        horizontalLayout_28->addWidget(IndepGroup);


        verticalLayout_7->addLayout(horizontalLayout_28);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        progressBar = new QProgressBar(centralwidget);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setValue(24);

        horizontalLayout_8->addWidget(progressBar);

        buttonSolve = new QPushButton(centralwidget);
        buttonSolve->setObjectName(QString::fromUtf8("buttonSolve"));

        horizontalLayout_8->addWidget(buttonSolve);


        verticalLayout_7->addLayout(horizontalLayout_8);


        horizontalLayout_29->addLayout(verticalLayout_7);

        SolverWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(SolverWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 826, 20));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menubar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        SolverWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(SolverWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        SolverWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addSeparator();
        menuFile->addAction(actionSave_Simulation);
        menuFile->addAction(actionSave_CSV);
        menuFile->addSeparator();
        menuFile->addAction(actionQuit);
        menuHelp->addAction(actionAbout);

        retranslateUi(SolverWindow);

        toolBox->setCurrentIndex(0);
        comboLeftBC->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(SolverWindow);
    } // setupUi

    void retranslateUi(QMainWindow *SolverWindow)
    {
        SolverWindow->setWindowTitle(QApplication::translate("SolverWindow", "Finite Difference Solver", 0, QApplication::UnicodeUTF8));
        actionQuit->setText(QApplication::translate("SolverWindow", "Quit", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("SolverWindow", "About", 0, QApplication::UnicodeUTF8));
        actionSave_Simulation->setText(QApplication::translate("SolverWindow", "Save As", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("SolverWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionSave_CSV->setText(QApplication::translate("SolverWindow", "Save CSV", 0, QApplication::UnicodeUTF8));
        label_37->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Domain Width</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_38->setText(QApplication::translate("SolverWindow", "Number of Nodes", 0, QApplication::UnicodeUTF8));
        label_39->setText(QApplication::translate("SolverWindow", "Width", 0, QApplication::UnicodeUTF8));
        spinWidth->setSuffix(QApplication::translate("SolverWindow", " m", 0, QApplication::UnicodeUTF8));
        label_40->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Time</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_41->setText(QApplication::translate("SolverWindow", "Step Size", 0, QApplication::UnicodeUTF8));
        spinDt->setSuffix(QApplication::translate("SolverWindow", " sec", 0, QApplication::UnicodeUTF8));
        label_42->setText(QApplication::translate("SolverWindow", "End Time", 0, QApplication::UnicodeUTF8));
        spintEnd->setSuffix(QApplication::translate("SolverWindow", " sec", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(DomainTab), QApplication::translate("SolverWindow", "Domain Settings", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("SolverWindow", "Mass Fractions", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("SolverWindow", "Protein", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("SolverWindow", "Fat", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("SolverWindow", "Carbohydrates", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("SolverWindow", "Fiber", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("SolverWindow", "Ash", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("SolverWindow", "Water", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("SolverWindow", "Ice", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(MatTab), QApplication::translate("SolverWindow", "Material Properties", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Initial Condition</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("SolverWindow", "Temperature", 0, QApplication::UnicodeUTF8));
        spinTInit->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Left Boundary Condition</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("SolverWindow", "Type", 0, QApplication::UnicodeUTF8));
        comboLeftBC->clear();
        comboLeftBC->insertItems(0, QStringList()
         << QApplication::translate("SolverWindow", "Insulation/Symmetry", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SolverWindow", "Convection", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SolverWindow", "Convection (Multiple Temperatures)", 0, QApplication::UnicodeUTF8)
        );
        lbch->setText(QApplication::translate("SolverWindow", "h", 0, QApplication::UnicodeUTF8));
        spinHLeft->setSuffix(QApplication::translate("SolverWindow", " W/(m^2 K)", 0, QApplication::UnicodeUTF8));
        lbctext->setText(QApplication::translate("SolverWindow", "External Temp", 0, QApplication::UnicodeUTF8));
        spinTExtLeft->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        lbctexthot->setText(QApplication::translate("SolverWindow", "External Temp (hot)", 0, QApplication::UnicodeUTF8));
        spinTExtHotLeft->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        lbctextcold->setText(QApplication::translate("SolverWindow", "External Temp (Cold)", 0, QApplication::UnicodeUTF8));
        spinTExtColdLeft->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        lbctheat->setText(QApplication::translate("SolverWindow", "Heating Time", 0, QApplication::UnicodeUTF8));
        spinTHeatLeft->setSuffix(QApplication::translate("SolverWindow", " sec", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Right Boundary Condition</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("SolverWindow", "Type", 0, QApplication::UnicodeUTF8));
        comboRightBC->clear();
        comboRightBC->insertItems(0, QStringList()
         << QApplication::translate("SolverWindow", "Insulation/Symmetry", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SolverWindow", "Convection", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("SolverWindow", "Convection (Multiple Temperatures)", 0, QApplication::UnicodeUTF8)
        );
        rbch->setText(QApplication::translate("SolverWindow", "h", 0, QApplication::UnicodeUTF8));
        spinHRight->setSuffix(QApplication::translate("SolverWindow", " W/(m^2 K)", 0, QApplication::UnicodeUTF8));
        rbctext->setText(QApplication::translate("SolverWindow", "External Temp", 0, QApplication::UnicodeUTF8));
        spinTExtRight->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        rbctexthot->setText(QApplication::translate("SolverWindow", "External Temp (hot)", 0, QApplication::UnicodeUTF8));
        spinTExtHotRight->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        rbctextcold->setText(QApplication::translate("SolverWindow", "External Temp (Cold)", 0, QApplication::UnicodeUTF8));
        spinTExtColdRight->setSuffix(QApplication::translate("SolverWindow", " K", 0, QApplication::UnicodeUTF8));
        rbctheat->setText(QApplication::translate("SolverWindow", "Heating Time", 0, QApplication::UnicodeUTF8));
        spinTHeatRight->setSuffix(QApplication::translate("SolverWindow", " sec", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(TTab), QApplication::translate("SolverWindow", "Heat Conduction", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Initial Condition</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("SolverWindow", "Concentration", 0, QApplication::UnicodeUTF8));
        spinC1Init->setSuffix(QApplication::translate("SolverWindow", " mol/L", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Reaction Kinetics</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("SolverWindow", "Prefactor", 0, QApplication::UnicodeUTF8));
        spinA1->setSuffix(QApplication::translate("SolverWindow", " sec^-1", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("SolverWindow", "Activation Energy", 0, QApplication::UnicodeUTF8));
        spinEa1->setSuffix(QApplication::translate("SolverWindow", " J/mol", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(C1Tab), QApplication::translate("SolverWindow", "Product Reaction Rate", 0, QApplication::UnicodeUTF8));
        label_33->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Initial Condition</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_35->setText(QApplication::translate("SolverWindow", "Concentration", 0, QApplication::UnicodeUTF8));
        spinC2Init->setPrefix(QString());
        spinC2Init->setSuffix(QApplication::translate("SolverWindow", " mol/L", 0, QApplication::UnicodeUTF8));
        label_32->setText(QApplication::translate("SolverWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Sans Serif'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">Reaction Kinetics</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        label_34->setText(QApplication::translate("SolverWindow", "Prefactor", 0, QApplication::UnicodeUTF8));
        spinA2->setSuffix(QApplication::translate("SolverWindow", " sec^-1", 0, QApplication::UnicodeUTF8));
        label_36->setText(QApplication::translate("SolverWindow", "Activation Energy", 0, QApplication::UnicodeUTF8));
        spinEa2->setSuffix(QApplication::translate("SolverWindow", " J/mol", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(C2Tab), QApplication::translate("SolverWindow", "Bacteria Reaction Rate", 0, QApplication::UnicodeUTF8));
        DepGroup->setTitle(QApplication::translate("SolverWindow", "Dependant Variables", 0, QApplication::UnicodeUTF8));
        checkTemp->setText(QApplication::translate("SolverWindow", "Temperature", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_ACCESSIBILITY
        checkk->setAccessibleDescription(QString());
#endif // QT_NO_ACCESSIBILITY
        checkk->setText(QApplication::translate("SolverWindow", "Thermal Diffusivity", 0, QApplication::UnicodeUTF8));
        checkProduct->setText(QApplication::translate("SolverWindow", "Product Concentration", 0, QApplication::UnicodeUTF8));
        checkBact->setText(QApplication::translate("SolverWindow", "Bacteria Concentration", 0, QApplication::UnicodeUTF8));
        IndepGroup->setTitle(QApplication::translate("SolverWindow", "Independant Variable", 0, QApplication::UnicodeUTF8));
        radioTime->setText(QApplication::translate("SolverWindow", "Time", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("SolverWindow", "(Time Index)", 0, QApplication::UnicodeUTF8));
        radioLength->setText(QApplication::translate("SolverWindow", "Length", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("SolverWindow", "(Node Number)", 0, QApplication::UnicodeUTF8));
        buttonPlot->setText(QApplication::translate("SolverWindow", "Plot Solution", 0, QApplication::UnicodeUTF8));
        buttonSolve->setText(QApplication::translate("SolverWindow", "Solve", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("SolverWindow", "File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("SolverWindow", "Help", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class SolverWindow: public Ui_SolverWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SOLVER_H
