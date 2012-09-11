/****************************************************************************
** Meta object code from reading C++ file 'solver.h'
**
** Created: Tue Sep 11 02:10:54 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../gui/solver.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'solver.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Solver[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
       8,    7,    7,    7, 0x0a,
      26,    7,    7,    7, 0x0a,
      43,    7,    7,    7, 0x0a,
      60,    7,    7,    7, 0x0a,
      68,    7,    7,    7, 0x0a,
      78,    7,    7,    7, 0x0a,
      93,    7,    7,    7, 0x0a,
     108,    7,    7,    7, 0x0a,
     124,    7,    7,    7, 0x0a,
     145,    7,    7,    7, 0x0a,
     161,    7,    7,    7, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Solver[] = {
    "Solver\0\0quitApplication()\0loadSimulation()\0"
    "saveSimulation()\0about()\0saveCSV()\0"
    "changeLBC(int)\0changeRBC(int)\0"
    "setMaxNode(int)\0setMaxTIndex(double)\0"
    "solveProblems()\0plotSolution()\0"
};

void Solver::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Solver *_t = static_cast<Solver *>(_o);
        switch (_id) {
        case 0: _t->quitApplication(); break;
        case 1: _t->loadSimulation(); break;
        case 2: _t->saveSimulation(); break;
        case 3: _t->about(); break;
        case 4: _t->saveCSV(); break;
        case 5: _t->changeLBC((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->changeRBC((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->setMaxNode((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: _t->setMaxTIndex((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 9: _t->solveProblems(); break;
        case 10: _t->plotSolution(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Solver::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Solver::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_Solver,
      qt_meta_data_Solver, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Solver::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Solver::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Solver::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Solver))
        return static_cast<void*>(const_cast< Solver*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int Solver::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
