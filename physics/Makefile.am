if CXX_GNU
    COM_CXXFLAGS = -Wall -W -pedantic -Wno-deprecated-declarations
else
if CXX_INTEL
    COM_CXXFLAGS = -wd1682 -Wall
endif
endif

bin_PROGRAMS = latan-2pt-fit latan-2pt-pseudosc-fit

latan_fit_2pt_SOURCES  = fit-2pt.cpp
latan_fit_2pt_CXXFLAGS = $(COM_CXXFLAGS)
latan_fit_2pt_LDFLAGS  = -L../lib/.libs -lLatAnalyze

latan_fit_phys_SOURCES  = fit-phys.cpp fit-phys-env.cpp
latan_fit_phys_CXXFLAGS = $(COM_CXXFLAGS)
latan_fit_phys_LDFLAGS  = -L../lib/.libs -lLatAnalyze

latan_2pt_pseudosc_fit_SOURCES  = 2pt-pseudosc-fit.cpp
latan_2pt_pseudosc_fit_CXXFLAGS = $(COM_CXXFLAGS)
latan_2pt_pseudosc_fit_LDFLAGS  = -L../lib/.libs -lLatAnalyze

ACLOCAL_AMFLAGS = -I .buildutils/m4
