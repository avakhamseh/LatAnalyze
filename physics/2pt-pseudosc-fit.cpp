
//clang++ -std=c++0x ava-2pt-fit-general.cpp -L/Users/Ava/local/lib -lLatAnalyze -lLatCore -I/Users/Ava/local/include -o ava-2pt-fit-general.x

#include <LatCore/OptParser.hpp>
#include <LatAnalyze/CompiledModel.hpp>
#include <LatAnalyze/Io.hpp>
#include <LatAnalyze/MatSample.hpp>
#include <LatAnalyze/Math.hpp>
#include <LatAnalyze/MinuitMinimizer.hpp>
#include <LatAnalyze/NloptMinimizer.hpp>
#include <LatAnalyze/Plot.hpp>
#include <LatAnalyze/XYSampleData.hpp>
#include <map>
                    //////////////////////////////  BOOT DATA FOR DECAY CONSTANTS!
using namespace std;
using namespace Latan;

enum CorrIndex {pp = 0, pa = 1, ap = 2, aa = 3};    
//enum {nCorr = 4}
enum {m = 0, zp = 1, za = 2};

int main(int argc, char *argv[])
{
  // parse arguments /////////////////////////////////////////////////////////
  OptParser            opt;
  bool                 parsed, doPlot, doHeatmap, doCorr, fold;
  string               model, outFileName, outFmt;
  Index                ti, tf, shift, nPar, thinning, nCorr;
//  Index                ti_pa, ti_ap, ti_aa, tf_pa, tf_ap, tf_aa
  double               svdTol;
  Minimizer::Verbosity verbosity;
  

  opt.addOption("" , "ti"       , OptParser::OptType::value  , false,
                "initial fit time");
  opt.addOption("" , "tf"       , OptParser::OptType::value  , false,
                "final fit time");
  // opt.addOption("" , "ti_pa"       , OptParser::OptType::value  , true,
  //               "initial fit time for pa channel");
  // opt.addOption("" , "ti_ap"       , OptParser::OptType::value  , true,
  //               "initial fit time for pa channel");
  // opt.addOption("" , "ti_pa"       , OptParser::OptType::value  , true,
  //               "initial fit time for aa channel");
  // opt.addOption("" , "fi_pa"       , OptParser::OptType::value  , true,
  //               "final fit time for pa channel");
  // opt.addOption("" , "fi_ap"       , OptParser::OptType::value  , true,
  //               "final fit time for pa channel");
  // opt.addOption("" , "fi_pa"       , OptParser::OptType::value  , true,
  //               "final fit time for aa channel");
  //opt.addOption("nCorr" , "nCorrelator" , OptParser::OptType::value  , true,
  //              "number of different correlators");
  opt.addOption("" , "pp" , OptParser::OptType::value  , true,   
                "psudoscalar-pseudoscalar correlator");
  opt.addOption("" , "pa" , OptParser::OptType::value  , true,
                "psudoscalar-axial correlator");
  opt.addOption("" , "ap" , OptParser::OptType::value  , true,
                "axial-pseudoscalar correlator");
  opt.addOption("" , "aa" , OptParser::OptType::value  , true,
                "axial-axial correlator");
  opt.addOption("t" , "thinning", OptParser::OptType::value  , true,
                "thinning of the time interval", "1");
  opt.addOption("s", "shift"    , OptParser::OptType::value  , true,
                "time variable shift", "0");
  opt.addOption("m", "model"    , OptParser::OptType::value  , true,
                "fit model (exp|exp2|exp3|cosh|cosh2|cosh3|<interpreter code>)", "cosh");
  opt.addOption("" , "nPar"     , OptParser::OptType::value  , true,
                "number of model parameters for custom models "
                "(-1 if irrelevant)", "-1");
  opt.addOption("" , "svd"      , OptParser::OptType::value  , true,
                "singular value elimination threshold", "0.");
  opt.addOption("v", "verbosity", OptParser::OptType::value  , true,
                "minimizer verbosity level (0|1|2)", "0");
  opt.addOption("o", "output", OptParser::OptType::value  , true,
                "output file", "");
  opt.addOption("" , "uncorr"   , OptParser::OptType::trigger, true,
                "only do the uncorrelated fit");
  opt.addOption("" , "fold"   , OptParser::OptType::trigger, true,
                "fold the correlator");
  opt.addOption("p", "plot"     , OptParser::OptType::trigger, true,
                "show the fit plot");
  opt.addOption("h", "heatmap"  , OptParser::OptType::trigger, true,
                "show the fit correlation heatmap");
  opt.addOption("", "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() != 0) or opt.gotOption("help"))  // IF I PUT SIZE = 0 THEN BAD ALLOC
  {
   // if (!parsed or !(opt.getArgs().size() <= 4) or opt.gotOption("help"))
  //{
    cerr << "usage: " << argv[0] << " <options> <correlator file>" << endl;
    cerr << endl << "Possible options:" << endl << opt << endl;

    return EXIT_FAILURE;
  }




  ti           = opt.optionValue<Index>("ti");
  tf           = opt.optionValue<Index>("tf");
  // ti_pa        = opt.optionValue<Index>("ti_pa");
  // tf_pa        = opt.optionValue<Index>("tf_pa");
  // ti_ap        = opt.optionValue<Index>("ti_ap");
  // tf_ap        = opt.optionValue<Index>("tf_ap");
  // ti_aa        = opt.optionValue<Index>("ti_aa");
  // tf_aa        = opt.optionValue<Index>("tf_aa");
  //nCorr        = opt.optionValue<Index>("nCorr");
  thinning     = opt.optionValue<Index>("t");
  shift        = opt.optionValue<Index>("s");
  model        = opt.optionValue("m");
  nPar         = opt.optionValue<Index>("nPar");
  svdTol       = opt.optionValue<double>("svd");
  outFileName  = opt.optionValue<string>("o");
  doCorr       = !opt.gotOption("uncorr");
  fold         = opt.gotOption("fold");
  doPlot       = opt.gotOption("p");
  doHeatmap    = opt.gotOption("h");
  switch (opt.optionValue<unsigned int>("v"))
  {
  case 0:
    verbosity = Minimizer::Verbosity::Silent;
    break;
  case 1:
    verbosity = Minimizer::Verbosity::Normal;
    break;
  case 2:
    verbosity = Minimizer::Verbosity::Debug;
    break;
  default:
    cerr << "error: wrong verbosity level" << endl;
    return EXIT_FAILURE;
  }

  //map<CorrIndex, string> corrFileName(nCorr);     
map<CorrIndex, string> corrFileName; //map does not get a size

  if (opt.gotOption("pp"))
  {
      corrFileName[pp] = opt.optionValue("pp");
  }
  if (opt.gotOption("pa"))
  {
      corrFileName[pa] = opt.optionValue("pa");
  }
  if (opt.gotOption("ap"))
  {
      corrFileName[ap] = opt.optionValue("ap");
  }
  if (opt.gotOption("aa"))
  {
      corrFileName[aa] = opt.optionValue("aa");
  }
  // corrType.size();
  // cout << "the size is " << corrType.size() << endl;

  //containers, iterators, range loop , struct

  // for (const auto &c: corrType)
  // {
  //    cout << c.first << ":" << c.second << '\n';
  // }

  //corrFileName = opt.getArgs().front();


  

  map<string, int> corrMap;
  corrMap["pp"] = 0;
  corrMap["pa"] = 1;
  corrMap["ap"] = 2;
  corrMap["aa"] = 3;

  for (const auto &c : corrMap)                      // & CHECK LOCATION!!!!!!!!!!
  {
    cout << c.first << ":" << c.second << endl;
  }

  // load correlator /////////////////////////////////////////////////////////
//     DMatSample tmp_pp, tmp_pa, tmp_ap, tmp_aa, corr_pp, corr_pa, corr_ap, corr_aa;
  Index      nSample, nt;

  map<CorrIndex, DMatSample> corr;
  DMatSample                 tmp;

  for (auto &c: corrFileName)
  {
    tmp     = Io::load<DMatSample>(c.second);
    nSample = tmp.size();
    nt      = tmp[central].rows();
    tmp     = tmp.block(0, 0, nt, 1);
    corr[c.first] = tmp;
  }


// for (Index i = 0, i < corrType.size(), i++)
// if (!corrType[i].empty()) { ..... }                   //IF IT'S VALUE IS NON ZERO THEN DO

   for (auto &c: corr)
  {
    tmp = corr[c.first];
    FOR_STAT_ARRAY(corr[c.first], s)
    {
      for (Index t = 0; t < nt; ++t)
      {
        corr[c.first][s]((t - shift + nt) % nt) = tmp[s](t);
      }
    }

    int sgn;
    if (fold)
    {
      if (c.first == aa || c.first == pp) {
        sgn = 1;
      }
      else {
        sgn = -1;
      }
      FOR_STAT_ARRAY(corr[c.first], s)
      {
        for (Index t = 0; t < nt; ++t)
        {
          corr[c.first][s](t) = 0.5 * (tmp[s](t) + sgn * tmp[s]((nt - t) % nt));

        }
      }

    }

  }

//

//     // make model //////////////////////////////////////////////////////////////
//     DoubleModel modCosh1, modCosh2, modCosh3, ...;
  DoubleModel modCosh_pp, modCosh_aa, modSinh_ap;
  vector<const DoubleModel *> modVec;                // WHY DID WE NOT DO A MAP HERE AS WELL

//     bool        coshModel = false;

  nPar = 3;

  modCosh_pp.setFunction([nt](const double * x, const double * p)     // discuss the instance .setFunction
  {
    return p[1] * p[1] * (exp(-p[0] * x[0]) + exp(-p[0] * (nt - x[0])));
  }, 1, nPar);

  modCosh_aa.setFunction([nt](const double * x, const double * p)
  {
    return p[2] * p[2] * (exp(-p[0] * x[0]) + exp(-p[0] * (nt - x[0])));
  }, 1, nPar);

  modSinh_ap.setFunction([nt](const double * x, const double * p)
  {
    return -p[1] * p[2] * (exp(-p[0] * x[0]) - exp(-p[0] * (nt - x[0])));
  }, 1, nPar);

  // modVec[pp] = &modCosh_pp;
  // modVec[ap] = &modSinh_ap;
  // modVec[pa] = &modSinh_ap;
  // modVec[aa] = &modCosh_aa;
  for (auto &c: corrFileName)
  {
    switch (c.first)
    {
      case pp:
        modVec.push_back(&modCosh_pp);
        break;
      case pa:
        modVec.push_back(&modSinh_ap);
        break;
      case ap:
        modVec.push_back(&modSinh_ap);
        break;
      case aa:
        modVec.push_back(&modCosh_aa);
        break;

    }
  }





//     if ((model == "exp") or (model == "exp1"))
//     {
//         nPar = 2;
//         mod.setFunction([](const double *x, const double *p)
//                         {
//                             return p[1]*exp(-p[0]*x[0]);
//                         }, 1, nPar);
//     }
//     else if (model == "exp2")
//     {
//         nPar = 4;
//         mod.setFunction([](const double *x, const double *p)
//                         {
//                             return p[1]*exp(-p[0]*x[0]) + p[3]*exp(-p[2]*x[0]);
//                         }, 1, nPar);
//     }
//     else if (model == "exp3")
//     {
//         nPar = 6;
//         mod.setFunction([](const double *x, const double *p)
//                         {
//                             return p[1]*exp(-p[0]*x[0]) + p[3]*exp(-p[2]*x[0])
//                                    + p[5]*exp(-p[4]*x[0]);
//                         }, 1, nPar);
//     }

//     else if ((model == "cosh") or (model == "cosh1"))
//     {
//         coshModel = true;
//         nPar      = 2;
//         modCosh1.setFunction([nt](const double *x, const double *p)
//                         {
//                             return p[1]*(exp(-p[0]*x[0])+exp(-p[0]*(nt-x[0])));
//                         }, 1, nPar);
//     }
//     else if (model == "cosh2")
//     {
//         coshModel = true;
//         nPar      = 4;
//         mod.setFunction([nt](const double *x, const double *p)
//                         {
//                             return p[1]*(exp(-p[0]*x[0])+exp(-p[0]*(nt-x[0])))
//                                  + p[3]*(exp(-p[2]*x[0])+exp(-p[2]*(nt-x[0])));
//                         }, 1, nPar);
//     }
//     else if (model == "cosh3")
//     {
//         coshModel = true;
//         nPar      = 6;
//         mod.setFunction([nt](const double *x, const double *p)
//                         {
//                             return p[1]*(exp(-p[0]*x[0])+exp(-p[0]*(nt-x[0])))
//                                  + p[3]*(exp(-p[2]*x[0])+exp(-p[2]*(nt-x[0])))
//                                  + p[5]*(exp(-p[2]*x[0])+exp(-p[4]*(nt-x[0])));
//                         }, 1, nPar);
//     }
//     else
//     {
//         if (nPar > 0)
//         {
//             mod = compile(model, 1, nPar);
//         }
//         else
//         {
//             cerr << "error: please specify the number of model parameter"
//                     " using the --nPar function" << endl;

//             return EXIT_FAILURE;
//         }
//     }

//     for (int c = 0; c < ???.size(); ++c)
//     {
//         if (???[c].type == ap)
//         {
//             modVec[i] = &modCosh1;
//         }
//         else
//         {
//             modVec[i] = &modSinh1;
//         }
//     }

//     // fit /////////////////////////////////////////////////////////////////////
  DMatSample          tvec(nSample);
  XYSampleData        data(nSample);
  SampleFitResult     fit;
  DVec                init(nPar);
  NloptMinimizer      globMin(NloptMinimizer::Algorithm::GN_CRS2_LM);
  MinuitMinimizer     locMin;
  vector<Minimizer *> unCorrMin{&globMin, &locMin};

  FOR_STAT_ARRAY(tvec, s)
  {
    tvec[s] = DVec::LinSpaced(nt, 0, nt - 1);
  }
  data.addXDim(nt, "t/a", true);

for (auto &c: corrFileName)
  {
    switch (c.first)
    {
      case pp:
        data.addYDim("C_pp(t)");
        break;
      case pa:
        data.addYDim("C_pa(t)");
        break;
      case ap:
        data.addYDim("C_ap(t)");
        break;
      case aa:
        data.addYDim("C_aa(t)" );
        break;

    }
  }

  //for (auto &c: corr)
  //{ 
  //if (!corr[c.first].empty()) { ..... } 

  data.setUnidimData(tvec, corr[pp], corr[pa], corr[ap], corr[aa]); //*****************************************//
  // DEBUG_VAR(corr[pp][central](0));

  modCosh_pp.parName().setName(m, "m");
  modCosh_pp.parName().setName(zp, "Z_p");
  modCosh_pp.parName().setName(za, "Z_a");

  init(m) = log(data.y(nt / 4, pp)[central] / data.y(nt / 4 + 1, 0)[central]);    // WHICH OF THE Y-DATA POINTS ???
  init(zp) = sqrt(data.y(nt / 4, pp)[central] / (exp(-init(0) * nt / 4)));        // data.y(const Index k, const Index j)  CHECK THIS!
  init(za) = sqrt(data.y(nt / 4, aa)[central] / (exp(-init(0) * nt / 4)));        // Row-major order, k = time, j = y-data but k can be generalise in the case of global

//  DEBUG_MAT(init);



  globMin.setLowLimit(m, 0.);
  globMin.setHighLimit(m, 10.*init(m));
  globMin.setLowLimit(zp, -10.*init(zp));
  globMin.setHighLimit(zp, 10.*init(zp));
  globMin.setLowLimit(za, -10.*init(za));
  globMin.setHighLimit(za, 10.*init(za));
  locMin.setLowLimit(m, 0.);


  globMin.setPrecision(0.001);
  globMin.setMaxIteration(100000);
  globMin.setVerbosity(verbosity);
  locMin.setMaxIteration(1000000);
  locMin.setVerbosity(verbosity);

  for (unsigned int i = 0; i < nCorr; ++i)
    for (Index t = 0; t < nt; ++t)
    {
      data.fitPoint((t >= ti) and (t <= tf)
                    and ((t - ti) % thinning == 0), t, i);
    }
  if (doCorr)
  {
    cout << "-- uncorrelated fit..." << endl;
  }
  //cout << "using model '" << model << "'" << endl;
  data.setSvdTolerance(svdTol);
  data.assumeYYCorrelated(false, 0, 0);
  fit = data.fit(unCorrMin, init, modVec);
  fit.print();
  if (doCorr)
  {
    cout << "-- correlated fit..." << endl;
    //cout << "using model '" << model << "'" << endl;
    init = fit[central];
    for (unsigned int i1 = 0; i1 < nCorr; ++i1)
      for (unsigned int i2 = 0; i2 < nCorr; ++i2)
      {
        data.assumeYYCorrelated(true, i1, i2);
      }
    fit = data.fit(locMin, init, modVec);
    fit.print();
  }

  // plots ///////////////////////////////////////////////////////////////////
  if (doPlot)
  {
    Plot       p;
   // vector<DMatSample> effMass(nCorr);
    map<enum CorrIndex, DMatSample> effMass;           
    DVec       effMassT, fitErr;
    Index      maxT =  (nt - 2);
    double     e0, e0Err;

    p << PlotRange(Axis::x, 0, nt - 1);
    p << LogScale(Axis::y);
    p << Color("rgb 'blue'") << PlotPredBand(fit.getModel(_, pp), 0, nt - 1);
    p << Color("rgb 'blue'") << Title("PP") << PlotFunction(fit.getModel(central, pp), 0, nt - 1);

    p << Color("rgb 'red'") << PlotPredBand(fit.getModel(_, pa), 0, nt - 1);
    p << Color("rgb 'red'") << Title("PA") << PlotFunction(fit.getModel(central, pa), 0, nt - 1);

    p << Color("rgb 'orange'") << PlotPredBand(fit.getModel(_, ap), 0, nt - 1);
    p << Color("rgb 'orange'") << Title("AP") << PlotFunction(fit.getModel(central, ap), 0, nt - 1);

    p << Color("rgb 'green'") << PlotPredBand(fit.getModel(_, aa), 0, nt - 1);
    p << Color("rgb 'green'") << Title("AA") << PlotFunction(fit.getModel(central, aa), 0, nt - 1);

    p << Color("rgb 'black'") << PlotData(data.getData(), 0, pp);
    p << Color("rgb 'black'") << PlotData(data.getData(), 0, pa);
    p << Color("rgb 'black'") << PlotData(data.getData(), 0, ap);
    p << Color("rgb 'black'") << PlotData(data.getData(), 0, aa);
    p.display();

for (auto &c: effMass)
  {
    effMass[c.first].resize(nSample);
    effMass[c.first].resizeMat(maxT, 1); 
  }


    effMassT.setLinSpaced(maxT, 1, maxT);
    fitErr = fit.variance().cwiseSqrt();
    e0     = fit[central](0);
    e0Err  = fitErr(0);
    //     if (coshModel)
    //     {
    
    for (auto &c: effMass)
  {
      FOR_STAT_ARRAY(effMass[c.first], s)
      {
        for (Index t = 1; t < nt - 1; ++t)
        {
          if (c.first == aa || c.first == pp) {
            effMass[c.first][s](t - 1) = acosh((corr[c.first][s](t - 1) + corr[c.first][s](t + 1))
                                         / (2.*corr[c.first][s](t)));
          }
          else {
            effMass[c.first][s](t - 1) = acosh((corr[c.first][s](t - 1) + corr[c.first][s](t + 1))
                                         / (2.*corr[c.first][s](t)));
          }
        }
      }
    }

    //     }

    //     else
    //     {
    //         FOR_STAT_ARRAY(effMass, s)
    //         {
    //             for (Index t = 1; t < nt; ++t)
    //             {
    //                 effMass[s](t - 1) = log(corr[s](t-1)/corr[s](t));
    //             }
    //         }
    //     }
    p.reset();
    p << PlotRange(Axis::x, 1, maxT);
    p << PlotRange(Axis::y, e0 - 40.*e0Err, e0 + 40.*e0Err);
    p << Color("rgb 'blue'") << PlotBand(0, maxT, e0 - e0Err, e0 + e0Err);
    p << Color("rgb 'blue'") << PlotHLine(e0);
    p << Color("rgb 'black'") << PlotData(effMassT, effMass[pp]);
    p << Color("rgb 'red'") << PlotData(effMassT, effMass[pa]);
    p << Color("rgb 'orange'") << PlotData(effMassT, effMass[ap]);
    p << Color("rgb 'green'") << PlotData(effMassT, effMass[aa]);
    
    //DEBUG_MAT(effMass[pp][central]);
   
    p.display();
    //p.save("test");
    // }
    // if (doHeatmap)
    // {
    //     Plot  p;
    //     Index n  = data.getFitVarMat().rows();
    //     DMat  id = DMat::Identity(n, n);

    //     p << PlotMatrix(Math::varToCorr(data.getFitVarMat()));
    //     p << Caption("correlation matrix");
    //     p.display();
    //     if (svdTol > 0.)
    //     {
    //         p.reset();
    //         p << PlotMatrix(id - data.getFitVarMat()*data.getFitVarMatPInv());
    //         p << Caption("singular space projector");
    //         p.display();
    //     }
  }

  // output //////////////////////////////////////////////////////////////////
  if (!outFileName.empty())
  {
    Io::save(fit, outFileName);
  }

  return EXIT_SUCCESS;
}

