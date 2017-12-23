
inline double sign(double a, double b)
{
  return (b >= 0.0 ? fabs(a) : -fabs(a));
}


inline double max(double a, double b)
{
  return (a >= b ? a : b);
}

inline double min(double a, double b)
{
  return (a <= b ? a : b);
}

inline int sign(int a, int b)
{
  return (b >= 0 ? abs(a) : -abs(a));
}

inline int max(int a, int b)
{
  return (a >= b ? a : b);
}

inline int min(int a, int b)
{
  return (a <= b ? a : b);
}

