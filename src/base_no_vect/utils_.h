double norm1est_(double *A, int n, int m);

double norm_1(double *A, int n);

void m_abs(double *A, double *result, int n);

void mcm(double *A, double *result, double scalar, int n);

void m_sum(double *A, double *B, double *result, int n);

void matrix_multiply_blas(int m, int n, int k, double *A, double *B, double **C);

int is_triangular(double *A, int n);

void solve_system(double *A, double *b, double **x, int n);

double *allocate_matrix(int n, int m);

double *identity_matrix(int n);

const double theta[14] = {0, 0, 0, 1.495585217958292e-2, 0, 2.539398330063230e-1,
                              0, 9.504178996162932e-1, 0, 2.097847961257068e0,
                              0, 0, 0, 4.25};

// Matrix b definition
const double b3[] = {120.0, 60.0, 12.0, 1.0};
const double b5[] = {30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};
const double b7[] = {17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0};
const double b9[] = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0};
const double b13[] = {
        64764752532480000.0,
        32382376266240000.0,
        7771770303897600.0,
        1187353796428800.0,
        129060195264000.0,
        10559470521600.0,
        670442572800.0,
        33522128640.0,
        1323241920.0,
        40840800.0,
        960960.0,
        16380.0,
        182.0,
        1.0
    };

const double c[] = {0, 0, 0, 100800.0, 0 , 10059033600.0, 0, 4487938430976000.0, 0, 5914384781877411840000.0, 0,0,0, 113250775606021113483283660800000000.0};

double *A2;
double *A4;
double *A6;
double *A8;
