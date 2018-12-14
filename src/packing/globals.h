#include <exception>
#include <string>
#include <vector>

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

class ClashRejection : public std::exception {
  const char* what() const throw() { return "Found an unreoslvable clash of shapes."; }
} ClashRejection;

#endif /* GLOBALS_H */
