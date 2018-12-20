#include <cmath>

#ifndef GLOBALS_H
#define GLOBALS_H

#define MAXVAR 30
#define MAXREPLICAS 500
#define MAXITER 20
#define MAXGRAPHICSIZE 1000
#define GRAPHICFRAME 100
#define CENTRE_SPOT 1
#define AXES_SHOWN 1

#define EPS 1e-9

#define ALLOWFLIPS 1
#define MAXSTEPS 360000
#define CYCLES 32
#define POLYGON_REP 1

#define WARNINGS 0

double
temperature_distribution(double old_val, double new_val, double kT, size_t replicas) {
  return std::exp(
      ((1.0 / old_val - 1.0 / new_val) / kT) + replicas * std::log(old_val / new_val));
}


#endif /* GLOBALS_H */
