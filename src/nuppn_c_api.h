#ifndef IMYANN_NUPPN_C_API_H
#define IMYANN_NUPPN_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

int nuppn_init(void);
int nuppn_num_species(void);
int nuppn_get_species(int index, char *name, int name_len, int *z, int *a);
int nuppn_integrate(double rho0, double temp0_k, double rho1, double temp1_k,
                    double *xnuc, int count, double dt, double *dedt);
int nuppn_integrate_to_time(double rho, double temp_k, double *xnuc, int count,
                            double final_time, double initial_dt,
                            double max_dt, double dt_factor, int max_steps,
                            double *dedt);
int nuppn_last_substeps(void);
int nuppn_history_size(void);
int nuppn_get_history_step(int index, double *time, double *rho,
                           double *temp_k, double *dt, double *dedt,
                           int *substeps, double *xnuc, int count);
int nuppn_compute_nse(double rho, double temp_k, double ye, double *xnuc,
                      int count);
int nuppn_get_dxdt(double rho, double temp_k, const double *xnuc, int count,
                   double *dxdt_out);
int nuppn_evaluate_reactions(double rho, double temp_k, const double *xnuc,
                             int count);
int nuppn_num_reactions(void);
int nuppn_get_reaction(int index, int *in1, int *in1_count, int *in2,
                       int *in2_count, int *out1, int *out1_count, int *out2,
                       int *out2_count, double *rate, double *flow,
                       double *q_value, int *weak);

#ifdef __cplusplus
}
#endif

#endif
