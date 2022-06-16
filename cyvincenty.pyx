import cython
# from cython.cimports.libc.math import sin, cos, tan, atan, atan2, pow, fabs, sqrt

# WGS 84
cdef int a = 6378137  # meters
cdef double f = 1 / 298.257223563
cdef double b = 6356752.314245  # meters; b = (1 - f)a
cdef double PI = 3.14159265358979323846264338327950288

cdef double MILES_PER_KILOMETER = 0.621371

cdef int MAX_ITERATIONS = 200
cdef double CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001

cdef extern from "math.h":
    double sin(double x)

cdef extern from "math.h":
    double cos(double x)

cdef extern from "math.h":
    double tan(double x)

cdef extern from "math.h":
    double atan(double x)

cdef extern from "math.h":
    double atan2(double x, double y)

cdef extern from "math.h":
    double pow(double x, double y)

cdef extern from "math.h":
    double fabs(double x)

cdef extern from "math.h":
    double sqrt(double x)

cdef double degreesToRadians(double angleDegrees):
    return ((angleDegrees) * PI / 180.0)

cdef class Distance:
    cdef double lat1, long1, lat2, long2, km, _km, miles, meters


    def __cinit__(self, double lat1, long1, lat2, long2) -> None:
        self.lat1 = lat1
        self.long1 = long1
        self.lat2 = lat2
        self.long2 = long2
        self._km = self.vincenty()

    @property
    def km(self):
        return self._km

    @property
    def miles(self):
        return self._km * MILES_PER_KILOMETER

    @property
    def meters(self):
        return self._km * 1000.0

    def vincenty(self):
        """
        Vincenty's formula (inverse method) to calculate the distance (in
        kilometers or miles) between two points on the surface of a spheroid
        Doctests:
        >>> vincenty((0.0, 0.0), (0.0, 0.0))  # coincident points
        0.0
        >>> vincenty((0.0, 0.0), (0.0, 1.0))
        111.319491
        >>> vincenty((0.0, 0.0), (1.0, 0.0))
        110.574389
        >>> vincenty((0.0, 0.0), (0.5, 179.5))  # slow convergence
        19936.288579
        >>> vincenty((0.0, 0.0), (0.5, 179.7))  # failure to converge
        >>> boston = (42.3541165, -71.0693514)
        >>> newyork = (40.7791472, -73.9680804)
        >>> vincenty(boston, newyork)
        298.396057
        >>> vincenty(boston, newyork, miles=True)
        185.414657
        """

        # short-circuit coincident points
        if self.lat1 == self.lat2 and self.long1 == self.long2:
            return 0.0

        cdef double U1 = atan((1 - f) * tan(degreesToRadians(self.lat1)))
        cdef double U2 = atan((1 - f) * tan(degreesToRadians(self.lat2)))
        cdef double L = degreesToRadians(self.long2 - self.long1)
        cdef double Lambda = L

        cdef double sinU1 = sin(U1)
        cdef double cosU1 = cos(U1)
        cdef double sinU2 = sin(U2)
        cdef double cosU2 = cos(U2)

        cdef double sinLambda, cosLambda, sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM, C, LambdaPrev

        for _ in range(MAX_ITERATIONS):
            sinLambda = sin(Lambda)
            cosLambda = cos(Lambda)
            sinSigma = sqrt(pow(cosU2 * sinLambda, 2.0) +
                                 pow(cosU1 * sinU2 - sinU1 * cosU2 * cosLambda, 2))
            if sinSigma == 0:
                return 0.0  # coincident points
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
            sigma = atan2(sinSigma, cosSigma)
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
            cosSqAlpha = 1 - pow(sinAlpha, 2)
            try:
                cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
            except ZeroDivisionError:
                cos2SigmaM = 0
            C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
            LambdaPrev = Lambda
            Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                                   (cos2SigmaM + C * cosSigma *
                                                    (-1 + 2 * pow(cos2SigmaM, 2))))
            if fabs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
                break  # successful convergence
        else:
            return None  # failure to converge

        cdef double uSq = cosSqAlpha * (pow(a, 2) - pow(b, 2)) / (pow(b, 2))
        cdef double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
        cdef double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
        cdef double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                     (-1 + 2 * pow(cos2SigmaM, 2)) - B / 6 * cos2SigmaM *
                     (-3 + 4 * pow(sinSigma, 2)) * (-3 + 4 * pow(cos2SigmaM, 2))))
        cdef double s = b * A * (sigma - deltaSigma)

        s /= 1000.0  # meters to kilometers
        return s
