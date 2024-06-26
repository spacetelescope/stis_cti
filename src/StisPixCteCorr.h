/* constants describing the CTE parameters reference file */
#define NUM_PHI 12  /* number of phi values in cte params file */
#define NUM_PSI 17  /* number of psi nodes in cte params file (also # of rows in table) */
#define NUM_LOGQ 4  /* number of log q columns in psi array */
#define NUM_LEV 298 /* number of specified Q levels */
#define NUM_SCALE 2 /* number of time dependant CTE scale points */

/* constants describing the CTE characterization */
#define MAX_TAIL_LEN 60 /* CTE trails are characterized out to 100 pixels */
#define MAX_PHI 150000    /* max number of phi nodes */
#define CTE_REF_ROW 1024  /* row from which CTE is measured */

/* constants for different instrument names */
#define ACSWFC 1

/* parameters of readout noise decomposition routines */
#define NOISE_MODEL 1

/* function prototypes */
double CalcCteFrac(const double expstart, const double scalemjd[NUM_SCALE],
                   const double scaleval[NUM_SCALE]);
int InterpolatePsi(const double chg_leak[NUM_PSI*NUM_LOGQ], const int psi_node[NUM_PSI],
                   double chg_leak_interp[MAX_TAIL_LEN*NUM_LOGQ],
                   double chg_open_interp[MAX_TAIL_LEN*NUM_LOGQ]);
int InterpolatePhi(const double dtde_l[NUM_PHI], const int q_dtde[NUM_PHI],
                   const int shft_nit, double dtde_q[MAX_PHI]);
int FillLevelArrays(const double chg_leak_kt[MAX_TAIL_LEN*NUM_LOGQ],
                    const double chg_open_kt[MAX_TAIL_LEN*NUM_LOGQ],
                    const double dtde_q[MAX_PHI], const int levels[NUM_LEV],
                    double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
                    double chg_open_lt[MAX_TAIL_LEN*NUM_LEV],
                    double dpde_l[NUM_LEV]);
int DecomposeRN(const int arrx, const int arry, const double data[ /* arrx*arry */],
                const double read_noise, const int noise_model,
                double sig_arr[ /* arrx*arry */ ], double noise_arr[ /* arrx*arry */ ]);
int FixYCte(const int arrx, const int arry, const double sig_cte[ /* arrx*arry */ ],
            double sig_cor[ /* arrx*arry */ ], const int sim_nit, const int shft_nit,
            const double too_low, double cte_frac[ /* arrx*arry */ ],
            const int levels[NUM_LEV], const double dpde_l[NUM_LEV],
            const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
            const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);
int AddYCte(const int arrx, const int arry, const double sig_cte[ /* arrx*arry */ ],
            double sig_cor[ /* arrx*arry */ ], const int shft_nit,
            const double cte_frac[ /* arrx*arry */ ], const int levels[NUM_LEV],
            const double dpde_l[NUM_LEV],
            const double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV],
            const double chg_open_lt[MAX_TAIL_LEN*NUM_LEV]);
