#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
struct particle { float x, y, z; };
const float G = 6.67e-11;

double wtime()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void calculate_forces_line(struct particle *p, struct particle *f, float *m, int n)
{
  for (int i = 0; i < n - 1; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      // Вычисление силы, действующей на тело i со стороны j
      float dist = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2) +
      powf(p[i].z - p[j].z, 2));
      float mag = (G * m[i] * m[j]) / powf(dist, 2);
      struct particle dir = {
        .x = p[j].x - p[i].x,
        .y = p[j].y - p[i].y,
        .z = p[j].z - p[i].z
      };
      // Сумма сил, действующих на тело i
      f[i].x += mag * dir.x / dist;
      f[i].y += mag * dir.y / dist;
      f[i].z += mag * dir.z / dist;
      // Сумма сил, действующих на тело j (симметричность)
      f[j].x -= mag * dir.x / dist;
      f[j].y -= mag * dir.y / dist;
      f[j].z -= mag * dir.z / dist;
    }
  }
}

void calculate_forces(struct particle *p, struct particle *f, float *m, int n)
{
  #pragma omp parallel for
  for (int i = 0; i < n - 1; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      // Вычисление силы, действующей на тело i со стороны j
      float dist = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2) +
      powf(p[i].z - p[j].z, 2));
      float mag = (G * m[i] * m[j]) / powf(dist, 2);
      struct particle dir = {
        .x = p[j].x - p[i].x,
        .y = p[j].y - p[i].y,
        .z = p[j].z - p[i].z
      };
      #pragma omp critical
      {
        // Сумма сил, действующих на тело i
        f[i].x += mag * dir.x / dist;
        f[i].y += mag * dir.y / dist;
        f[i].z += mag * dir.z / dist;
        // Сумма сил, действующих на тело j (симметричность)
        f[j].x -= mag * dir.x / dist;
        f[j].y -= mag * dir.y / dist;
        f[j].z -= mag * dir.z / dist;
      }
    }
  }
}

void move_particles(struct particle *p, struct particle *f, struct particle *v, float *m, int n,
double dt)
{
  for (int i = 0; i < n; i++)
  {
    struct particle dv = {
      .x = f[i].x / m[i] * dt,
      .y = f[i].y / m[i] * dt,
      .z = f[i].z / m[i] * dt,
    };
    struct particle dp = {
      .x = (v[i].x + dv.x / 2) * dt,
      .y = (v[i].y + dv.y / 2) * dt,
      .z = (v[i].z + dv.z / 2) * dt,
    };
    v[i].x += dv.x;
    v[i].y += dv.y;
    v[i].z += dv.z;
    p[i].x += dp.x;
    p[i].y += dp.y;
    p[i].z += dp.z;
    f[i].x = f[i].y = f[i].z = 0;
  }
}

void func_line(struct particle *p, struct particle *f, struct particle *v, float *m, int n, double dt, double *tforces, double *tmove)
{
  for (double t = 0; t <= 1; t += dt) // Цикл по времени (модельному)
  {
    *tforces -= wtime();
    calculate_forces_line(p, f, m, n); // Вычисление сил – O(N^2)
    *tforces += wtime();
    *tmove -= wtime();
    move_particles(p, f, v, m, n, dt); // Перемещение тел O(N)
    *tmove += wtime();
  }
}

void func(struct particle *p, struct particle *f, struct particle *v, float *m, int n, double dt, double *tforces, double *tmove, int threads)
{
  for (double t = 0; t <= 1; t += dt) // Цикл по времени (модельному)
  {
    *tforces -= wtime();
    #pragma omp parallel num_threads(threads)
    {
      calculate_forces(p, f, m, n); // Вычисление сил – O(N^2)
    }
    *tforces += wtime();
    *tmove -= wtime();
    move_particles(p, f, v, m, n, dt); // Перемещение тел O(N)
    *tmove += wtime();
  }
}

void init(struct particle *p, struct particle *f, struct particle *v, float *m, int n)
{
  for (int i = 0; i < n; i++)
  {
    p[i].x = rand() / (float)RAND_MAX - 0.5;
    p[i].y = rand() / (float)RAND_MAX - 0.5;
    p[i].z = rand() / (float)RAND_MAX - 0.5;
    v[i].x = rand() / (float)RAND_MAX - 0.5;
    v[i].y = rand() / (float)RAND_MAX - 0.5;
    v[i].z = rand() / (float)RAND_MAX - 0.5;
    m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
    f[i].x = f[i].y = f[i].z = 0;
  }
}

void copy(struct particle *p, struct particle *f, struct particle *v, float *m, struct particle *p_cp, struct particle *f_cp, struct particle *v_cp, float *m_cp, int n)
{
  for(int i = 0; i < n; i++)
  {
    p[i].x = p_cp[i].x;
    p[i].y = p_cp[i].y;
    p[i].z = p_cp[i].z;
    v[i].x = v_cp[i].x;
    v[i].y = v_cp[i].y;
    v[i].z = v_cp[i].z;
    m[i] = m_cp[i];
    f[i].x = f_cp[i].x;
    f[i].y = f_cp[i].y;
    f[i].z = f_cp[i].z;
  }
}

int main(int argc, char *argv[])
{
  FILE *fp = fopen("critical.dat", "w");
  double boost = 0;
  double ttotal, tinit = 0, tforces = 0, tmove = 0, t_line = 0, t_parallel = 0;
  ttotal = wtime();
  int n = (argc > 1) ? atoi(argv[1]) : 10;
  char *filename = (argc > 2) ? argv[2] : NULL;
  char buf[10] = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '\0'}, bufn[5] = { ' ', ' ', ' ', ' ', '\0'};
  tinit = -wtime();
  struct particle *p = malloc(sizeof(*p) * n); // Положение частиц (x, y, z)
  struct particle *f = malloc(sizeof(*f) * n); // Сила, действующая на каждую частицу (x, y, z)
  struct particle *v = malloc(sizeof(*v) * n); // Скорость частицы (x, y, z)
  float *m = malloc(sizeof(*m) * n); // Масса частицы
  struct particle *p_cp = malloc(sizeof(*p) * n);
  struct particle *f_cp = malloc(sizeof(*f) * n);
  struct particle *v_cp = malloc(sizeof(*v) * n);
  float *m_cp = malloc(sizeof(*m) * n);
  init(p, f, v, m, n);
  tinit += wtime();
  double dt = 1e-5;
  printf("# NBody (n=%d)\n", n);
  printf("# Init time (sec): %.6f\n", tinit);
  copy(p_cp, f_cp, v_cp, m_cp, p, f, v, m, n);
  t_line = wtime();
  func_line(p, f, v, m, n, dt, &tforces, &tmove);
  t_line = wtime() - t_line;
  printf("# Elapsed time (sec) for serial: %.6f\n", t_line);
  copy(p, f, v, m, p_cp, f_cp, v_cp, m_cp, n);
  for(int i = 2; i <= 8; i += 2)
  {
    func(p, f, v, m, n, dt, &tforces, &tmove, i);
    t_parallel = tforces + tmove;
    printf("# Elapsed time (sec) for %d threads: ttotal %.6f, tforces %.6f, tmove %.6f\n",
    i, t_parallel, tforces, tmove);
    if(i != 8)
      copy(p, f, v, m, p_cp, f_cp, v_cp, m_cp, n);
    boost = t_line / t_parallel;
    bufn[0] = i + '0';
    sprintf(buf, "%.6f", boost);
    buf[8] = '\n';
    fputs(bufn, fp);
    fputs(buf, fp);
  }
  ttotal = wtime() - ttotal;
  if (filename)
  {
    FILE *fout = fopen(filename, "w");
    if (!fout)
    {
      fprintf(stderr, "Can't save file\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++)
      fprintf(fout, "%15f %15f %15f\n", p[i].x, p[i].y, p[i].z);
    fclose(fout);
  }
  fclose(fp);
  free(p);
  free(f);
  free(v);
  free(m);
  free(p_cp);
  free(f_cp);
  free(v_cp);
  free(m_cp);
  return 0;
}