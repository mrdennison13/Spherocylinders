#include <stdio.h>
#include <gsl/gsl_rng.h>

gsl_rng * ran;
const gsl_rng_type * T;


void random2_num_int_()
{

  gsl_rng_env_setup();
  ran = gsl_rng_alloc (gsl_rng_taus);
  unsigned int seed;

  FILE *fp=fopen("/dev/urandom","r");
  fread(&seed,1,sizeof(unsigned int),fp);
  fclose(fp);
  gsl_rng_set(ran,seed);

}


double random2_num_()
{
  return gsl_rng_uniform (ran);
}




