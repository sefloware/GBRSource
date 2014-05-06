TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.cpp

HEADERS += \
    lincs/type.h \
    lincs/geometricConstraint.h \
    lincs/diffusion.h \
    stoRebounding.h \
    lincs/geometricConstraintBlas.h \
    lincs/lincsBlas.h \
    lincs/blasRef.h \
    lincs/lincs.h \
    bendingForce.h \
    obstacles.h \
    SCnanochannel.h \
    GCsphere.h \
    GCrod.h \
    GCplain.h \
    GCcylinder.h \
    SCsphere.h \
    common.h \
    Fanglespring1.h \
    Fanglespring2.h \
    rebound_PB.h \
    lincs_blas/type.h \
    lincs_blas/lincs.h \
    lincs_blas/geometricConstraint.h \
    lincs_blas/diffusion.h

QMAKE_CXXFLAGS += -O2

INCLUDEPATH += C:\Qt\Qt5.1.1\Tools\mingw48_32\include
